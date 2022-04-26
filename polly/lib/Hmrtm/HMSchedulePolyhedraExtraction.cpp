#include "polly/HMSchedulePolyhedraExtraction.h"

#define DEBUG_TYPE "hm-polyhedra-extraction"
/**
 * This implements the schedule space extraction.
 */

using namespace polly::hmrtm;


polly::hmrtm::StmtCoeffInfo::StmtCoeffInfo(const int itStart,
                                           const int nrIt,
                                           const int parStart,
                                           const int cstIdx,
                                           const isl::id &identifier) : itStart(itStart),
                                                                        nrIt(nrIt),
                                                                        parStart(parStart),
                                                                        cstIdx(cstIdx),
                                                                        identifier(
                                                                            identifier) {}

polly::hmrtm::DomainCoeffInfo::DomainCoeffInfo(const isl::union_set &domain_)
    : nrIt(compNrIt(domain_)), nrParPS(compNrParPS(domain_)),
      stmtInfo(compStmtInfo(domain_)), universe(compUniverse(domain_)),
      domain(domain_), nrStmts(compNrStmts(domain_)),
      dim(nrIt + nrParPS * nrStmts + nrStmts) {

}
int polly::hmrtm::DomainCoeffInfo::compNrIt(const isl::union_set &domain) {
  auto setList = domain.get_set_list();
  int nrIt = 0;
  for (int i = 0; i < setList.size(); ++i) {
    nrIt += setList.get_at(i).dim(isl::dim::set);
  }
  return nrIt;
}
int polly::hmrtm::DomainCoeffInfo::compNrParPS(const isl::union_set &domain) {
  return domain.params().dim(isl::dim::param);
}
int polly::hmrtm::DomainCoeffInfo::compNrStmts(const isl::union_set &domain) {
  return domain.n_set();
}
llvm::StringMap<polly::hmrtm::StmtCoeffInfo> polly::hmrtm::DomainCoeffInfo::compStmtInfo(
    const isl::union_set &domain) {
  llvm::SmallVector<isl::set, 0> stmtSets;
  auto setList = domain.get_set_list();
  for (int i = 0; i < setList.size(); ++i) {
    stmtSets.push_back(setList.get_at(i));
  }
  std::sort(stmtSets.begin(), stmtSets.end(), [](isl::set &a, isl::set &b) {
    return a.get_tuple_name() < b.get_tuple_name();
  });
  llvm::StringMap<StmtCoeffInfo> stmtInfo;
  int domainParDim = compNrParPS(domain);
  int nrStmts = compNrStmts(domain);
  int nrIt = compNrIt(domain);
  int count = 0;
  for (size_t i = 0; i < stmtSets.size(); ++i) {
    auto &stmt = stmtSets[i];
    int stmtNrIt = stmt.dim(isl::dim::set);
    auto result = stmtInfo.try_emplace(stmt.get_tuple_name(),
                                       count,
                                       stmtNrIt,
                                       nrIt + domainParDim * i,
                                       nrIt + domainParDim * nrStmts + i,
                                       stmt.get_tuple_id());
    if (!result.second) {
      LLVM_DEBUG({
                   llvm::dbgs() << "Cleary, something went wrong....";
                 });
    }
    count += stmtNrIt;
  }
  return stmtInfo;
}
isl::set polly::hmrtm::DomainCoeffInfo::compUniverse(const isl::union_set &domain) {
  int nrIt = compNrIt(domain);
  int nrStmts = compNrStmts(domain);
  int domainParDim = compNrParPS(domain);
  auto ctx = domain.get_ctx();
  auto space = isl::manage(isl_space_set_alloc(ctx.get(), 0,
                                               static_cast<unsigned int>(
                                                   nrIt + nrStmts
                                                       * (domainParDim
                                                           + 1))));
  return isl::set::universe(space);
}
polly::hmrtm::Dependence::Dependence(const isl::basic_map &baseMap_,
                                     const isl::basic_set &weakConstr_,
                                     const isl::basic_set &strongConstr_) :
    baseMap(baseMap_), weakConstr(weakConstr_), strongConstr(strongConstr_) {
}
bool polly::hmrtm::Dependence::operator<(const polly::hmrtm::Dependence &other) const {
  return baseMap.to_str() < other.baseMap.to_str();
}

namespace {
/**
 * Separate the dependence polyhedron in individual dependencies.
 * @param deps A dependence polyhedron.
 * @param domain The associated domain.
 * @return A list of individual dependencies. When united, should be equal
 * to deps
 */
static llvm::SmallVector<isl::basic_map, 0> preprocess(isl::union_map deps,
                                                       isl::union_set domain) {

  llvm::SmallVector<isl::basic_map, 0> result;
  int maxSplit = 100;
  auto mapList = deps.get_map_list();

  // First split the union map into individual maps
  for (int j = 0, size = mapList.size(); j < size; ++j) {
    auto lMap = mapList.get_at(j);
    llvm::SmallVector<isl::basic_map, 10> intermediate;
    int i = 0;
    auto mapUnwrapped = lMap.range_is_wrapping()
                        ? lMap.uncurry().domain().unwrap()
                        : lMap;
    auto depsMap = mapUnwrapped;
    // Honestly, I do not exaclty why the lexmin operations are used, but
    // the code below is correct, as either the current map is split into
    // its basic parts by using lexmin and substracting it...
    do {
      auto dep = depsMap.lexmin().coalesce();
      auto basicMapList = dep.get_basic_map_list();
      for (int k = 0, sizeK = basicMapList.size(); k < sizeK; ++k) {
        auto basicMap = basicMapList.get_at(k);
        intermediate.push_back(basicMap);
        i = i + 1;
        if (i > maxSplit) {
          break;
        }
      }
      depsMap = depsMap.subtract(dep);
    } while (!depsMap.is_empty() && i <= maxSplit);
    if (i <= maxSplit) {
      result.append(intermediate.begin(), intermediate.end());
    } else {
      // ... or by simply adding all basic maps directly
      // Either way, the map just gets split into individual dependencies.
      auto basicMapList = mapUnwrapped.get_basic_map_list();
      for (int k = 0, sizeK = basicMapList.size(); k < sizeK; ++k) {
        result.push_back(basicMapList.get_at(k).remove_divs());
      }
    }
  }

  // remove empty constraints
  result.erase(std::remove_if(result.begin(), result.end(),
                              [&domain](const isl::basic_map &el) {
                                return isl::union_map(el).intersect_domain(
                                    domain).intersect_range(domain).is_empty();
                              }), result.end());

  // Sort by tuple names.
  std::sort(result.begin(), result.end(),
            [](isl::basic_map &first, isl::basic_map &second) {
              return std::max(first.get_tuple_name(isl::dim::in),
                              first.get_tuple_name(isl::dim::out)) <
                  std::max(second.get_tuple_name(isl::dim::in),
                           second.get_tuple_name(isl::dim::out));
            });

  return result;
}

/**
 * This computes a schedule constraint for a given dependency.
 * @param dep The dependency map.
 * @param domInfo The domain coefficient info.
 * @param strongSatisfy Whether or not to strongly satisfy the constraint.
 * @return The set describing the schedule constraint.
 */
static isl::basic_set compSchedConstrForDep(isl::basic_map dep,
                                            const DomainCoeffInfo &domInfo,
                                            bool strongSatisfy) {
  auto ctx = dep.get_ctx();
  isl::space oldSpace = dep.get_space();
  isl::id depSrcId = oldSpace.get_tuple_id(isl::dim::in);
  isl::id depDstId = oldSpace.get_tuple_id(isl::dim::out);
  // Check whether the dependency points to the same statement or not
  // This makes a difference as in the former case, only the schedule
  // coefficients of that particular statement are constrained, whereas
  // in the latter case, two statements are influenced by this.
  bool inEqOut = depSrcId.get_name() == depDstId.get_name();

  auto srcDim = dep.dim(isl::dim::in);
  auto dstDim = dep.dim(isl::dim::out);
  auto parDim = dep.dim(isl::dim::param);
  auto depPoly =
      dep.move_dims(isl::dim::in, srcDim, isl::dim::out, 0, dstDim).domain();
  // Coefficients applies the Farkas lemma to the dependency
  auto schedCoefficients = depPoly.coefficients().unwrap();

  auto njuDim = inEqOut ? srcDim : srcDim + dstDim;

  // it is short for iterator
  // the iterator variables of the src dimension have to be negated
  auto postprocIt =
      isl::map::universe(isl::space(ctx, 0, srcDim + dstDim, njuDim));
  for (decltype(srcDim) i = 0; i < srcDim; ++i) {
    postprocIt = postprocIt.oppose(isl::dim::in, i, isl::dim::out, i);
  }
  // if statements are not the same, leave dst dimensions untouched
  auto njuOff = inEqOut ? 0 : srcDim;
  for (decltype(dstDim) i = 0; i < dstDim; ++i) {
    postprocIt =
        postprocIt.equate(isl::dim::in, srcDim + i, isl::dim::out, njuOff + i);
  }
  schedCoefficients = schedCoefficients.apply_range(postprocIt.affine_hull());

//  LLVM_DEBUG({
//    llvm::dbgs() << "\nFor constraint: " << dep.to_str() << "\n";
//    llvm::dbgs() << "Sched coefficients in the middle: " << schedCoefficients.to_str() << "\n";
//  });

  isl::basic_map postprocPar(nullptr);
  // for same statement, they cancel out each other, so just introduce the positivity
  // constraint here
  if (inEqOut) {
    postprocPar = isl::basic_map::universe(isl::space(ctx, 0, parDim + 1, 0));
    auto zeroVal = isl::val::zero(ctx);
    auto valC = strongSatisfy ? isl::val::negone(ctx) : zeroVal;
    postprocPar = postprocPar.fix_val(isl::dim::in, 0, valC);
    for (decltype(parDim) i = 0; i < parDim; ++i) {
      postprocPar = postprocPar.fix_val(isl::dim::in, i + 1, zeroVal);
    }
    // otherwise, make sure that the difference of both is either >= 0 for weak
    // and >= 1 for strong
  } else {
    auto mAff = isl::multi_aff::zero(
        isl::space(ctx, parDim + 1, parDim + 1, 2 * parDim + 2));
    auto lSp = isl::local_space(isl::space(ctx, parDim + 1, parDim + 1));
    auto oneVal = isl::val::one(ctx);
    auto affC1 = isl::aff::var_on_domain(lSp, isl::dim::param, 0);
    auto affC2 = isl::aff::var_on_domain(lSp, isl::dim::set, 0).add(affC1);
    if (strongSatisfy) {
      affC2 = isl::aff(lSp, oneVal).add(affC2);
    }
    mAff = mAff.set_aff(0, affC1);
    mAff = mAff.set_aff(1, affC2);
    LLVM_DEBUG({
      llvm::dbgs() << "\nFor dep: " << dep.to_str() << "\n";
      llvm::dbgs() << "mAff: " << mAff.to_str() << "\n";
    });
    for (decltype(parDim) i = 0; i < parDim; ++i) {
      auto affP1 = isl::aff::var_on_domain(lSp, isl::dim::param, i + 1);
      auto
          affP2 = isl::aff::var_on_domain(lSp, isl::dim::set, i + 1).add(affP1);
      mAff = mAff.set_aff(i + 2, affP1);
      mAff = mAff.set_aff(i + 2 + parDim, affP2);
    }
    postprocPar =
        isl::basic_map::from_multi_aff(mAff).project_out(isl::dim::param,
                                                         0,
                                                         parDim + 1);
  }

  schedCoefficients = schedCoefficients.apply_domain(postprocPar);

//  LLVM_DEBUG({
//               llvm::dbgs() << "\nFor constraint: " << dep.to_str() << "\n";
//               llvm::dbgs() << "Postproc par: " << postprocPar.to_str() << "\n";
//               llvm::dbgs() << "Sched coefficients after apply_domain: " << schedCoefficients.to_str() << "\n";
//  });

  // These entries have to exist, so do not validate.
  // If the program crashes, we have a major bug anyways.
  auto stmtSrcInfo = domInfo.stmtInfo.find(depSrcId.get_name())->second;
  auto stmtDstInfo = domInfo.stmtInfo.find(depDstId.get_name())->second;
  auto srcItOff = stmtSrcInfo.itStart;
  auto srcNrIt = stmtSrcInfo.nrIt;
  auto srcParOff = stmtSrcInfo.parStart;
  auto srcCstIdx = stmtSrcInfo.cstIdx;
  auto dstItOff = stmtDstInfo.itStart;
  auto dstNrIt = stmtDstInfo.nrIt;
  auto dstParOff = stmtDstInfo.parStart;
  auto dstCstIdx = stmtDstInfo.cstIdx;

  auto nrPar = domInfo.nrParPS;
  auto solSpMap = isl::basic_map::from_domain_and_range(
      schedCoefficients.reverse().wrap().flatten(),
      isl::basic_set::universe(isl::space(ctx, 0, domInfo.dim)));
  std::remove_const<decltype(domInfo.dim)>::type off = 0;
  for (decltype(srcNrIt) i = 0; i < srcNrIt; ++i) {
    solSpMap =
        solSpMap.equate(isl::dim::in, off + i, isl::dim::out, srcItOff + i);
  }
  if (!inEqOut) {
    off += srcNrIt;
    for (decltype(dstNrIt) i = 0; i < dstNrIt; ++i) {
      solSpMap =
          solSpMap.equate(isl::dim::in, off + i, isl::dim::out, dstItOff + i);
    }
    off = off + dstNrIt;

    solSpMap = solSpMap.equate(isl::dim::in, off, isl::dim::out, srcCstIdx);
    ++off;
    solSpMap = solSpMap.equate(isl::dim::in, off, isl::dim::out, dstCstIdx);
    ++off;

    for (decltype(nrPar) i = 0; i < nrPar; ++i) {
      solSpMap =
          solSpMap.equate(isl::dim::in, off + i, isl::dim::out, srcParOff + i);
    }
    off += nrPar;
    for (decltype(nrPar) i = 0; i < nrPar; ++i) {
      solSpMap =
          solSpMap.equate(isl::dim::in, off + i, isl::dim::out, dstParOff + i);
    }
  }

  return solSpMap.range();
}

/**
 * Calculate schedule constraints for the dependencies.
 * @param bMaps
 * @param domInfo
 * @param computeOut isl timeout to use for this.
 * @return
 */
static llvm::SmallVector<Dependence,
                         0> bMaps2Deps(const llvm::SmallVector<isl::basic_map,
                                                               0> &bMaps,
                                       const DomainCoeffInfo &domInfo,
                                       unsigned long computeOut) {
  llvm::SmallVector<Dependence, 0> dependences;
  for (auto &bMap: bMaps) {
    auto ctx = bMap.get_ctx();
    auto oldMaxOps = isl_ctx_get_max_operations(ctx.get());
    isl_ctx_set_max_operations(ctx.get(), computeOut);
    auto dSimpl = bMap.remove_redundancies();
    // insert here
    isl::basic_set weakConstr = compSchedConstrForDep(dSimpl, domInfo, false);
    isl::basic_set strongConstr = compSchedConstrForDep(dSimpl, domInfo, true);

    isl_ctx_reset_operations(ctx.get());
    isl_ctx_set_max_operations(ctx.get(), oldMaxOps);
    dependences.emplace_back(bMap, weakConstr, strongConstr);
  }
  return dependences;
}
}

static void printFlow(const isl::union_flow& flow) {
  LLVM_DEBUG({
    llvm::dbgs() << "Must Dep: " << flow.get_must_dependence().to_str() << "\n";
    llvm::dbgs() << "May Dep: " << flow.get_may_dependence().to_str() << "\n";
    llvm::dbgs() << "Must No Src: " << flow.get_must_no_source().to_str() << "\n";
    llvm::dbgs() << "May No Src: " << flow.get_may_no_source().to_str() << "\n";
  });
}

std::pair<DomainCoeffInfo,
          llvm::SmallVector<Dependence,
                            0>> polly::hmrtm::calcDepsAndDomInfo(const isl::set &params,
                                                                 isl::union_set domain,
                                                                 isl::union_map sched,
                                                                 const isl::union_map &reads,
                                                                 const isl::union_map &writes) {

  // first, initialize an empty map.
  auto empty = isl::union_map::empty(domain.get_space().params());
  // intersect schedule with the domain to obtain the correct mapping.
  auto schedule = sched.intersect_domain(domain).coalesce();

  // Use the isl approximate flow analysis (described in
  // Presburger Formulas and Polyhedral Compilation by Verdoolaege, p.140ff.)
  // to calculate memory dependencies.
  // WAW and reads in between.
  isl::union_access_info writeFlow(writes);
  writes.dump();
  reads.dump();
  schedule.dump();
  writeFlow = writeFlow.set_must_source(writes);
  writeFlow = writeFlow.set_may_source(reads);
  writeFlow = writeFlow.set_schedule_map(schedule);
  auto flow = writeFlow.compute_flow();
  LLVM_DEBUG({
    llvm::dbgs() << "-------\nWrite flow analysis:\n";
    printFlow(flow);
  });

  // This also contains the must dependencies.
  auto antiOut = flow.get_may_dependence();

  // RAW
  isl::union_access_info readFlow(reads);
  readFlow = readFlow.set_must_source(writes);
  readFlow = readFlow.set_may_source(empty);
  readFlow = readFlow.set_schedule_map(schedule);
  flow = readFlow.compute_flow();
  LLVM_DEBUG({
    llvm::dbgs() << "-------\nRead flow analysis:\n";
    printFlow(flow);
  });
  // Unite both dependence polyhedra.
  auto deps = antiOut.unite(flow.get_must_dependence()).coalesce();
  LLVM_DEBUG(llvm::dbgs() << "Deps after flow analysis:\n" << deps.to_str() << "\n");
  DomainCoeffInfo domInfo(domain);
  // cannot use SmallVectorImpl as RVO can only work if the types are exactly
  // equal! This makes sense as the object used in the function is the one
  // defined here, but this can only work if the types match exactly.

  // The deps currently is a union map, but we need to separate those. This is
  // done in preprocess.
  llvm::SmallVector<isl::basic_map, 0> depList = preprocess(deps, domain);


  // Compute schedule space constraints from dependencies.
  llvm::SmallVector<Dependence, 0> dependences
      = bMaps2Deps(depList, domInfo, 0);

  return std::make_pair(std::move(domInfo), std::move(dependences));

}

llvm::SmallVector<isl::set, 0> polly::hmrtm::constructMultiDimensionalSchedulePolyhedrons(
    const isl::set &params,
    isl::union_set domain,
    isl::union_map sched,
    const isl::union_map &reads,
    const isl::union_map &writes,
    double prob2Carry,
    std::pair<DomainCoeffInfo,
              llvm::SmallVector<Dependence,
                                0>> &depsAndDomInfo) {
  DomainCoeffInfo& domInfo = depsAndDomInfo.first;
  llvm::SmallVector<Dependence, 0>& depsOld = depsAndDomInfo.second;
  llvm::SmallVector<Dependence, 0> deps;
  // remove duplicate dependencies
  for(auto &dep: depsOld) {
    bool toBeInserted = true;
    for(auto &depCmp: deps) {
      if(dep.strongConstr.to_str() == depCmp.strongConstr.to_str()) {
        toBeInserted = false;
        break;
      }
    }
    if(toBeInserted) {
      deps.push_back(dep);
    }
  }


  isl::set schedulePoly = domInfo.universe;
  for(auto &dep: deps) {
    schedulePoly = schedulePoly.intersect(dep.weakConstr);
  }
  if(schedulePoly.is_empty()) {
    LLVM_DEBUG(llvm::dbgs() << "Cannot satisfy the constraints!");
  }
  llvm::SmallVector<Dependence*, 0> depsOrder;
  llvm::SmallPtrSet<Dependence*, 8> uncarriedDeps;
  llvm::SmallPtrSet<Dependence*, 8> carriedInCurrentDim;
  llvm::SmallVector<Dependence*, 0> stronglyCarried;
  llvm::SmallVector<isl::set, 0> schedulePolys;
  for(auto &dep: deps) {
    depsOrder.push_back(&dep);
    uncarriedDeps.insert(&dep);
  }
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> uniformDouble(0.0, 1.0);
  std::shuffle(depsOrder.begin(), depsOrder.end(), rng);
  while(!depsOrder.empty()) {
    stronglyCarried.clear();
    if(uniformDouble(rng) <= prob2Carry) {
      auto depsLeft = depsOrder.size();

      // inclusive upper bound
      auto numToCarry = std::uniform_int_distribution<decltype(depsLeft)>(0, depsLeft)(rng);
      // TODO: remove later
      if(numToCarry < 1) {
        numToCarry = 1;
      }
      LLVM_DEBUG(llvm::dbgs() << "numToCarry: " << numToCarry << "\n");
      stronglyCarried.assign(depsOrder.begin(), depsOrder.begin() + numToCarry);
    }
    auto currentScheduleDim = domInfo.universe;
    // first, weakly carry all dependences
    for(Dependence* dep: uncarriedDeps) {
      currentScheduleDim = currentScheduleDim.intersect(dep->weakConstr);
    }
    // then, try to strongly carry selected dependences
    carriedInCurrentDim.clear();
    for(Dependence* dep: stronglyCarried) {
      auto currentScheduleDimMaybe
          = currentScheduleDim.intersect(dep->strongConstr);
      if(!currentScheduleDimMaybe.is_empty()) {
        LLVM_DEBUG(llvm::dbgs() << "Carrying " << dep->strongConstr.to_str()
                                << " strongly!\n");
        carriedInCurrentDim.insert(dep);
        currentScheduleDim = currentScheduleDimMaybe;
      }
    }
    LLVM_DEBUG(llvm::dbgs() << "Currently carried: "
                            << carriedInCurrentDim.size() << "\n");
    for(Dependence *dep: carriedInCurrentDim) {
      uncarriedDeps.erase(dep);
    }
    depsOrder.erase(std::remove_if(depsOrder.begin(),
                                   depsOrder.end(),
                                   [&carriedInCurrentDim] (Dependence * dep) {
                                     // this is the usual c++ set contains check
                                     // set::contains gets added in c++ 20
                                     return carriedInCurrentDim.find(dep) != carriedInCurrentDim.end();
                                   }), depsOrder.end());
    schedulePolys.push_back(currentScheduleDim);
  }
  return schedulePolys;
}
llvm::SmallVector<isl::set, 0> polly::hmrtm::constructMultiDimensionalSchedulePolyhedrons(
    const isl::set &params,
    isl::union_set domain,
    isl::union_map sched,
    const isl::union_map &reads,
    const isl::union_map &writes,
    double prob2Carry) {
  auto depsAndDomInfo = calcDepsAndDomInfo(params, domain, sched, reads, writes);
  return constructMultiDimensionalSchedulePolyhedrons(params, domain, sched, reads, writes, prob2Carry, depsAndDomInfo);
}
llvm::raw_ostream &llvm::operator<<(llvm::raw_ostream &out,
                                    const polly::hmrtm::StmtCoeffInfo &info) {
  return out << "StmtCoeffInfo"
             << "\n  itStart = " << info.itStart
             << "\n  nrIt = " << info.nrIt
             << "\n  parStart = " << info.parStart
             << "\n  cstIdx = " << info.cstIdx
             << "\n";
}
llvm::raw_ostream &llvm::operator<<(llvm::raw_ostream &out,
                                    const polly::hmrtm::DomainCoeffInfo &info) {
  out << "DomainCoeffInfo:"
      << "\nnrIt = " << info.nrIt
      << "\nnrParPS = " << info.nrParPS
      << "\nuniverse = " << info.universe.to_str()
      << "\ndomain = " << info.domain.to_str();
  for (auto &stmt: info.stmtInfo) {
    out << "\n" << stmt.first().str() << ": " << stmt.second;
  }
  return out << "\n";
}
void polly::hmrtm::testCalcDepsAndDomInfo(isl::ctx ctx) {
  isl::set params(ctx, "[n] -> { : 0 < n}");

  isl::union_set domain(ctx, "[n] -> { S[i,j] : 0 <= i < n and 0 <= j < n }");
  domain = domain.add_set(
      isl::set(ctx, "[n] -> { T[i, j] : 0 <=i < 2*n and 0 <= j < 2*n }"));

  isl::union_map schedule(ctx, "[n] -> { S[i, j] -> [i, j] }");
  schedule = schedule.add_map(isl::map(ctx, "[n] -> { T[i, j] -> [i, j] }"));

  isl::union_map writes(ctx, "[n] -> { S[i,j] -> A[i, j] }");
  writes = writes.add_map(isl::map(ctx, "[n] -> { T[i,j] -> B[i, j] }"));

  isl::union_map reads(ctx, "[n] -> { T[i,j] -> B[i, j - 1]}");
  reads = reads.add_map(isl::map(ctx, "[n] -> { S[i,j] -> A[i - 1, j] }"));

  DomainCoeffInfo domainCoeffInfo(domain);
  LLVM_DEBUG(llvm::dbgs() << "In test calc: " << domainCoeffInfo << "\n");
  //calcDepsAndDomInfo(params, domain, schedule, reads, writes);
  auto result = constructMultiDimensionalSchedulePolyhedrons(
      params, domain, schedule, reads, writes, 0.9);
  LLVM_DEBUG(llvm::dbgs() << "Schedule polyhedron: \n");
  for(auto &set: result) {
    LLVM_DEBUG(llvm::dbgs() << set.to_str() << "\n");
  }
}
