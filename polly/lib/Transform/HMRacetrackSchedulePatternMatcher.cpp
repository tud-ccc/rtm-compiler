//===------ HMRacetrackScheduleOptimizer.h ----------------------*- C++ -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file contains the implementation of the alternation transformation
// for the MMA. It first runs the isl scheduler,
// then uses the already existing recognition code for MMA and applies the
// schedule change when the MMA is detected.
//
//===----------------------------------------------------------------------===//
#include "polly/HMScheduleShiftCalculator.h"
#include "polly/HMScheduleTreeTikzPrinter.h"
#include "polly/DependenceInfo.h"
#include "polly/LinkAllPasses.h"
#include "polly/ScheduleOptimizer.h"
#include "polly/ScopInfo.h"
#include "polly/ScheduleTreeTransform.h"
#include "polly/Support/GICHelper.h"

#include "llvm/Support/CommandLine.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ADT/StringMap.h"

#include "isl/ctx.h"
#include "isl/options.h"
#include "isl/printer.h"
#include "isl/schedule.h"
#include "isl/schedule_node.h"
#include "isl/union_map.h"
#include "isl/union_set.h"

#define DEBUG_TYPE "polly-opt-hm-pattern"


using namespace llvm;
using namespace polly;

static cl::opt<bool>
    HMPollyEnablePatternMatching("polly-hm-match-pattern", cl::desc("Enable the pattern matching."),
                            cl::init(true), cl::ZeroOrMore);

static cl::opt<bool>
    HMPollyEnablePlutoBeforePatterns("polly-hm-pluto-before-pattern", cl::desc("Enable Pluto before the pattern matching."),
                                 cl::init(true), cl::ZeroOrMore);

static cl::opt<bool>
    HMPollyEnableFullPermutationSearch("polly-hm-full-permutation-search", cl::desc("Try all loop permutations. Requires Pluto to be enabled."),
                                     cl::init(false), cl::ZeroOrMore);

static cl::opt<bool>
    HMPollyUseUnifiedPatternMatcher("polly-hm-unified-pattern-matcher", cl::desc("Enable the unified pattern matcher instead of the single ones."),
        cl::init(true), cl::ZeroOrMore);

static cl::opt<int> HMPollyPlutoOptions("polly-hm-pluto-options", cl::desc("Potential pluto options."), cl::init(2 | 32 | 64 | 128 | 256 | 512));

static cl::opt<std::string>
    HMPollyTraceDir("polly-hm-trace-dir",
              cl::desc("The directory for the trace files."),
              cl::value_desc("Directory path"), cl::ValueRequired,
              cl::init(""));

static cl::opt<std::string>
    HMPollyKernelName("polly-hm-kernel-name",
                    cl::desc("The name of this kernel."),
                    cl::ValueRequired,
                    cl::init(""));

namespace {
class HMRacetrackSchedulePatternMatcher : public ScopPass {
public:
  static char ID;

  explicit HMRacetrackSchedulePatternMatcher() : ScopPass(ID) {}

  /// Optimize the schedule of the SCoP @p S for SPM utilizing racetrack memory.
  bool runOnScop(Scop &S) override;

  /// Register all analyses and transformations required
  void getAnalysisUsage(AnalysisUsage &AU) const override;


};
} // namespace

// Basically, I need the walk that is done by the ScheduleOptimizer. It only
// allows to optimize the inner-most band nodes.
/**
 * This uses isl_schedule_node_map_descendant_bottom_up to visit and transform
 * every node in the given schedule. It provides the dependencies as an
 * additional argument to each call.
 * @param D The dependences
 * @param Schedule The schedule to transform.
 * @param fn The transformation function for each node.
 * @return The modified schedule.
 */
static isl::schedule walkSchedule(const Dependences &D, isl::schedule &Schedule,
    __isl_give isl_schedule_node *(*fn)(
        __isl_take isl_schedule_node *, void *)) {
  auto Root = Schedule.get_root();
  Root = isl::manage(isl_schedule_node_map_descendant_bottom_up(
      Root.release(), fn, const_cast<void *>(static_cast<const void *>(&D))
  ));
  return Root.get_schedule();
}

/**
 * Prints the C code for the current schedule of a scop. For each statement,
 * a statement function is called.
 * @param S The scop to print.
 */
static void printAstForScop(Scop &S) {
  LLVM_DEBUG({
    auto ast = isl::ast_build::from_context(S.getContext());
    auto printer = isl_printer_to_str(ast.get_ctx().get());
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
    auto scheduleTree = S.getScheduleTree();
    auto astNode = ast.node_from_schedule(scheduleTree);
    printer = isl_printer_print_ast_node(printer, astNode.get());
    std::string out(isl_printer_get_str(printer));
    dbgs() << out << "\n";
    isl_printer_free(printer);
  });
}

// Copied from the polly schedule optimizer pass. I need this.
/// Permute two dimensions of the band node.
///
/// Permute FirstDim and SecondDim dimensions of the Node.
///
/// @param Node The band node to be modified.
/// @param FirstDim The first dimension to be permuted.
/// @param SecondDim The second dimension to be permuted.
static isl::schedule_node permuteBandNodeDimensions(isl::schedule_node Node,
                                                    unsigned FirstDim,
                                                    unsigned SecondDim) {
  assert(isl_schedule_node_get_type(Node.get()) == isl_schedule_node_band &&
      isl_schedule_node_band_n_member(Node.get()) >
          std::max(FirstDim, SecondDim));
  auto PartialSchedule =
      isl::manage(isl_schedule_node_band_get_partial_schedule(Node.get()));
  auto PartialScheduleFirstDim = PartialSchedule.get_union_pw_aff(FirstDim);
  auto PartialScheduleSecondDim = PartialSchedule.get_union_pw_aff(SecondDim);
  PartialSchedule =
      PartialSchedule.set_union_pw_aff(SecondDim, PartialScheduleFirstDim);
  PartialSchedule =
      PartialSchedule.set_union_pw_aff(FirstDim, PartialScheduleSecondDim);
  Node = isl::manage(isl_schedule_node_delete(Node.release()));
  return Node.insert_partial_schedule(PartialSchedule);
}


static isl::map filteredMultiUnionPwAffToMap(
    const isl::multi_union_pw_aff &Aff, isl::id TupleId) {
  // TODO: this should merge multiple to a single map
  auto UnionAff = isl::union_pw_multi_aff(Aff);
  auto AffList = UnionAff.get_pw_multi_aff_list();
  isl::map ResultMap;
  for(int i = 0; i < AffList.size(); ++i) {
    auto SingleAff = AffList.get_at(i);
    if(SingleAff.get_tuple_id(isl::dim::in).to_str() == TupleId.to_str()) {
      ResultMap = isl::map::from_pw_multi_aff(SingleAff);
    } else {
      LLVM_DEBUG(dbgs() << "Hedging!\n");
    }
  }
  return ResultMap;
}

static isl::map createAlternatingPartialSchedule(
    const isl::map &Base, const isl::map &Alternating,
    int Sign, int Direction) {
  auto Schedule = Base.flat_range_product(Alternating); // S[i, j] -> [i], S[i, j] -> [j] => S[i, j] -> [i, j]
  Schedule = Schedule.add_dims(isl::dim::out, 2); // S[i, j] -> [i, j, o1, o2]
  auto Space = isl::local_space(Schedule.get_space());
  // modulus constraint
  auto Constraint = isl::constraint::alloc_equality(Space); // 1 i - 2 o1 = 0
  Constraint = Constraint.set_coefficient_si(isl::dim::out, 0, 1);
  Constraint = Constraint.set_coefficient_si(isl::dim::out, 2, -2);
  Constraint = Constraint.set_constant_si(Sign);
  Schedule = Schedule.add_constraint(Constraint);
  // direction constraint
  Constraint = isl::constraint::alloc_equality(Space); // j - o2 = 0
  Constraint = Constraint.set_coefficient_si(isl::dim::out, 1, -1);
  Constraint = Constraint.set_coefficient_si(isl::dim::out, 3, Direction);
  Schedule = Schedule.add_constraint(Constraint);

  // project out to get existential quantifier and to remove additional output
  Schedule = Schedule.project_out(isl::dim::out, 0, 3);
  LLVM_DEBUG(dbgs() << "New schedule: " << Schedule.to_str() << "\n");
  return Schedule;
}

static isl_bool nodeGetCoincident(const isl::schedule_node &Node, int Pos) {
  return isl_schedule_node_band_member_get_coincident(Node.get(), Pos);
}

static isl_bool nodeGetPermutable(const isl::schedule_node &Node) {
  return isl_schedule_node_band_get_permutable(Node.get());
}

static isl::union_set nodeGetAstOpts(const isl::schedule_node &Node) {
  return isl::manage(isl_schedule_node_band_get_ast_build_options(Node.get()));
}

static isl::schedule_node nodeSetPermutable(isl::schedule_node Node,
    isl_bool Permutable) {
  return isl::manage(
      isl_schedule_node_band_set_permutable(Node.release(), Permutable));
}

static int nodeBandMembers(isl::schedule_node Node) {
  return isl_schedule_node_band_n_member(Node.get());
}

static bool isNodeType(
    const isl::schedule_node &Node, isl_schedule_node_type Type) {
  return isl_schedule_node_get_type(Node.get()) == Type;
}

static isl::multi_union_pw_aff bandGetPartialSchedule(
    const isl::schedule_node &Node) {
  // TODO: what about get_partial_schedule_union_map?
  return isl::manage(isl_schedule_node_band_get_partial_schedule(Node.get()));
}


static int getSingleNonZeroIndex(const isl::set &Set) {
  auto Zero = isl::val::zero(Set.get_ctx());
  auto OptimizeCounter = 0;
  auto LastIndex = -1;
  for(int i = 0; i < Set.n_dim(); ++i) {
    auto Val = Set.plain_get_val_if_fixed(isl::dim::set, i);
    if(!Val.eq(Zero)) {
      ++OptimizeCounter;
      LastIndex = i;
    }
  }
  if(OptimizeCounter > 1) {
    LastIndex = -1;
  }
  return LastIndex;
}

static bool involveSameDimensions(const isl::map &First, const isl::map &Second,
    const isl::dim Dim) {
  if(!First.get_space().is_equal(Second.get_space())) {
    return false;
  }
  for(unsigned i = 0; i < First.dim(Dim); ++i) {
    if(First.involves_dims(Dim, i, 1) != Second.involves_dims(Dim, i, 1)) {
      return false;
    }
  }
  return true;
}

static isl::boolean isCoincidentLoop(const isl::union_map &Deps,
                                     const isl::union_map &Schedule) {
  auto DeltaSet = Deps.apply_domain(Schedule).apply_range(Schedule)
      .deltas().detect_equalities().coalesce();
  return DeltaSet.is_empty() || DeltaSet.is_equal(isl::set(Deps.get_ctx(), "{[0]}"));

}

static void groupMemoryAccessesByName(const ScopStmt * Stmt,
    std::vector<isl::map> &MemoryAccessGroups) {
  using MemAccessCollection = std::vector<MemoryAccess *>;

  // Group mem accesses to same reference together
  llvm::StringMap<MemAccessCollection> GroupedMemAccesses;
  for(auto *MemAccess: *Stmt) {
    auto StringId = MemAccess->getAccessRelation().get_tuple_id(isl::dim::out).to_str();
    if(GroupedMemAccesses.find(StringId) == GroupedMemAccesses.end()) {
      GroupedMemAccesses.insert(
          std::pair<StringRef, MemAccessCollection>(
              StringId, MemAccessCollection()));
    }
    GroupedMemAccesses.find(StringId)->second.push_back(MemAccess);
  }
  for(auto &MemAccessTuple: GroupedMemAccesses) {
    auto &MemAccesses = MemAccessTuple.second;
    auto MemAccessIt = MemAccesses.begin();
    auto Collected = (*MemAccessIt)->getAccessRelation();
    auto FirstRelation = Collected;
    bool IsEqualLoopOrder = true;
    ++MemAccessIt;
    // This does two things: First, it unites all the accesses into a single
    // map. Second, it checks whether all memory accesses have the same
    // base underlying access structure, i.e. [i + x][j+x]..., and not
    // something like A[i][j] and A[j][k], because in this case, we do
    // not know yet how to optimize.
    // Of course, if there is only one memory access, the pattern is
    // given, as it can only be validated by two structurally different
    // accesses.
    while(MemAccessIt != MemAccesses.end()) {
      auto CompareRelation = (*MemAccessIt)->getAccessRelation();
      Collected = Collected.unite(CompareRelation);
      auto Deltas = FirstRelation.apply_domain(CompareRelation).deltas();
      IsEqualLoopOrder = IsEqualLoopOrder && Deltas.is_singleton();
      ++MemAccessIt;
    }
    // only process if this condition is satisfied.
    if(IsEqualLoopOrder) {
      MemoryAccessGroups.emplace_back(Collected);
    }
  }
}

static void groupMemoryAccessesByNameAndEqualLoopOrder(
    const ScopStmt * Stmt,
    std::vector<std::vector<isl::map>> &MemoryAccessGroups) {
  using MemAccessCollection = std::vector<MemoryAccess *>;
  // Group mem accesses to same reference together
  llvm::StringMap<MemAccessCollection> GroupedMemAccesses;
  for(auto *MemAccess: *Stmt) {
    auto StringId = MemAccess->getAccessRelation().get_tuple_id(isl::dim::out).to_str();
    if(GroupedMemAccesses.find(StringId) == GroupedMemAccesses.end()) {
      GroupedMemAccesses.insert(
          std::pair<StringRef, MemAccessCollection>(
              StringId, MemAccessCollection()));
    }
    GroupedMemAccesses.find(StringId)->second.push_back(MemAccess);
  }
  for(auto &MemAccessTuple: GroupedMemAccesses) {
    auto &MemAccesses = MemAccessTuple.second;
    std::vector<std::vector<isl::map>> SameAccessOrder;
    std::vector<isl::map> UnitedMaps;
    for(auto &MemAccess: MemAccesses) {
      bool FoundEqualLoopOrder = false;
      auto CurrentRelation = MemAccess->getAccessRelation();
      for(auto i = 0u; i < SameAccessOrder.size(); ++i) {
        auto FirstRelation = SameAccessOrder[i][0];
        auto Deltas = FirstRelation.apply_domain(CurrentRelation).deltas();
        if(Deltas.is_singleton()) {
          FoundEqualLoopOrder = true;
          SameAccessOrder[i].push_back(CurrentRelation);
        }
      }
      if(!FoundEqualLoopOrder) {
        SameAccessOrder.emplace_back(std::vector<isl::map>{CurrentRelation});
      }
    }
    for(auto i = 0u; i < SameAccessOrder.size(); ++i) {
      for(auto j = 1u; j < SameAccessOrder[i].size(); ++j) {
        SameAccessOrder[i][0] = SameAccessOrder[i][0].unite(SameAccessOrder[i][j]);
      }
      UnitedMaps.push_back(SameAccessOrder[i][0]);
    }
    MemoryAccessGroups.emplace_back(std::move(UnitedMaps));
  }
}

static isl::set createAlternationBaseAstOption(const isl::map &AlternationBase) {
  auto Shifted = AlternationBase.move_dims(
      isl::dim::out, 0,
      isl::dim::in, 0,
      AlternationBase.dim(isl::dim::in));
  return Shifted.range().set_tuple_id(AlternationBase.get_tuple_id(isl::dim::in));
}

/**
 * Splits a band node consisting of multiple bands at the last one.
 * @param Node The band node to split. Calling this on any other node has no effect.
 * @param User Unused user pointer as second argument that allows to use this in
 * isl_schedule_node_map_descendant_bottom_up
 * @return The node for the last band.
 */
static __isl_give isl_schedule_node * splitBandNodes(
    __isl_take isl_schedule_node *Node, void *User) {
  if(isl_schedule_node_get_type(Node) == isl_schedule_node_band) {
    int BandMembers = isl_schedule_node_band_n_member(Node);
    if(BandMembers > 1) {
      Node = isl_schedule_node_band_split(Node, BandMembers - 1);
      isl_schedule_node * Child = isl_schedule_node_get_child(Node, 0);

      isl_schedule_node_free(Node);
      return Child;
    }
  }
  return Node;
}

static __isl_give isl_schedule_node * markAsCoincident(
    __isl_take isl_schedule_node *Node, void *User
    ) {
  if(isl_schedule_node_get_type(Node) == isl_schedule_node_band) {
    auto D = static_cast<Dependences *>(User);
    int ValidityKinds = Dependences::TYPE_RAW | Dependences::TYPE_WAR | Dependences::TYPE_WAW;
    isl::union_map Validity = D->getDependences(ValidityKinds);
    auto ManagedNode = isl::manage(Node);
    auto AffSched = bandGetPartialSchedule(ManagedNode);
    auto Sched = isl::union_map::from(AffSched);
    auto Coincident = isCoincidentLoop(Validity, Sched);
    Node = isl_schedule_node_band_member_set_coincident(
        ManagedNode.release(),
        0, Coincident.is_true() ? 1 : 0);
  }
  return Node;
}


static __isl_give isl_schedule_node * optimizeUnifiedPattern(
    __isl_take isl_schedule_node *Node, void *User) {
  // only operate on leafs again
  if(isl_schedule_node_get_type(Node) == isl_schedule_node_leaf) {
    LLVM_DEBUG(dbgs() << "Unified in leaf!\n");
    auto ManagedNode = isl::manage(Node);

    auto D = static_cast<Dependences *>(User);
    int ValidityKinds = Dependences::TYPE_RAW | Dependences::TYPE_WAR | Dependences::TYPE_WAW;
    isl::union_map Validity = D->getDependences(ValidityKinds);

    // current node is leaf, so at least transition to parent
    auto FilterNode = ManagedNode.parent();
    // find filter from which to descend later
    while(!isNodeType(FilterNode, isl_schedule_node_filter)
        && !isNodeType(FilterNode, isl_schedule_node_domain)) {
      FilterNode = FilterNode.parent();
    }
    // now find statement
    auto SubtreeScheduleUnion = FilterNode.get_subtree_schedule_union_map();
    // Make sure the subtree really only contains a single statement
    // TODO: is this really necessary? not here any more
    if(SubtreeScheduleUnion.n_map() == 1) {
      auto Domains = ManagedNode.get_domain();
      auto SubtreeSchedule = SubtreeScheduleUnion.get_map_list().get_at(0);
      auto StmtId = SubtreeSchedule.get_tuple_id(isl::dim::in);
      LLVM_DEBUG(dbgs() << "In stmt: " << StmtId.to_str() << "\n");
      // The stmt id stores the ScopStmt object
      auto *Stmt = static_cast<const ScopStmt *>(StmtId.get_user());
      using MemAccessCollection = std::vector<MemoryAccess *>;

      std::vector<std::vector<isl::map>> MemoryAccessGroups;
      groupMemoryAccessesByNameAndEqualLoopOrder(Stmt, MemoryAccessGroups);

      for(auto &Group: MemoryAccessGroups) {
        auto United = Group[0];
        for(auto i = 1u; i < Group.size(); ++i) {
          United = United.unite(Group[i]);
        }

        if (United.dim(isl::dim::out) > 0 && !United.is_injective()) {
          for (auto &Collected: Group) {
            std::vector<int> UninvolvedIndices;
            for (int i = 0; i < Collected.dim(isl::dim::in); ++i) {
              bool InvolvesDim = Collected.involves_dims(isl::dim::in, i, 1);
              if (!InvolvesDim) {
                Collected = Collected.fix_si(isl::dim::in, i, 0);
                UninvolvedIndices.push_back(i);
              }
            }

            // Now calculate the access distance set, i.e. that contains all
            // distances to the first access
            // Especially, the first access itself is represented as [0,0,...,0],

            auto Domain = Domains.extract_set(Collected.domain().get_space());
            LLVM_DEBUG({
                         dbgs() << "\nDomain: " << Domain.to_str()
                                << "\nCollection: " << Collected.to_str() << "\n";
                       });
            LLVM_DEBUG({
                         //dbgs() << "Special lexmin: " << Collected.reverse().lexmin().to_str() << "\n";
                       });
            auto Reversed = Collected.reverse();
            auto RawCtx = Reversed.get_ctx().get();
            // TODO: revert error checking afterwards
            auto IslOnError = isl_options_get_on_error(RawCtx);
            isl_options_set_on_error(RawCtx, ISL_ON_ERROR_CONTINUE);
            isl_ctx_reset_error(RawCtx);
            auto AccessLexMin = Reversed.lexmin();
            auto LexMinError = isl_ctx_last_error(RawCtx);
            if (LexMinError) {
              LLVM_DEBUG(dbgs() << "Got lexmin error: " << isl_ctx_last_error_msg(RawCtx) << "\n");
              LLVM_DEBUG(
                  dbgs() << "This is not critical, but it disables the pattern matching for this memory access!\n");
              continue;
            }
            for (auto Index: UninvolvedIndices) {
              Collected = Collected.drop_constraints_involving_dims(isl::dim::in, Index, 1);
              // Well, as the they are not involved, applying wont change anything...
              //AccessLexMin = AccessLexMin.drop_constraints_involving_dims(isl::dim::out, Index, 1);
            }
            Reversed = Collected.reverse();

            // TODO: add fail check for lexmin
            auto DeltaSet = Reversed.apply_domain(AccessLexMin).deltas();

            //std::vector<int> MultipleAccessLoops;
            // Next, create a schedule function for all of those accesses.
            auto BaseScheduleSpace = isl::space(DeltaSet.get_ctx(),
                                                DeltaSet.dim(isl::dim::param), DeltaSet.n_dim(), 1);
            BaseScheduleSpace = BaseScheduleSpace.set_tuple_id(isl::dim::in, DeltaSet.get_tuple_id());
            auto BaseSchedule = isl::map::universe(BaseScheduleSpace);

            // TODO: what if no dimension is involved?
            // One could possibly move to the next index, e.g. A[i][j][0]
            // With the default layout, this will provide no benefit, but with others, this might
            // find the access that is fastest changing, e.g. j in A[i][j]
            auto OnlyInnerLoopAccessMap = Reversed.drop_constraints_involving_dims(
                isl::dim::in, 0, Reversed.dim(isl::dim::in) - 1);
            for (int i = 0; i < OnlyInnerLoopAccessMap.dim(isl::dim::out); ++i) {
              bool InvolvesDim = OnlyInnerLoopAccessMap.involves_dims(isl::dim::out, i, 1);
              if (!InvolvesDim) {
                OnlyInnerLoopAccessMap = OnlyInnerLoopAccessMap.fix_si(isl::dim::out, i, 0);
              }
            }
            OnlyInnerLoopAccessMap = OnlyInnerLoopAccessMap.drop_constraints_involving_dims(
                isl::dim::in, Reversed.dim(isl::dim::in) - 1, 1);
            auto OnlyInnerLoopAccess = OnlyInnerLoopAccessMap.range();

            auto InnermostAccessIndex = getSingleNonZeroIndex(OnlyInnerLoopAccess);
            if (InnermostAccessIndex != -1) {
              auto AlternatedLoop = BaseSchedule.equate(
                  isl::dim::in, InnermostAccessIndex, isl::dim::out, 0);
              std::vector<isl::map> MultipleAccessSchedules;
              auto Dims = DeltaSet.n_dim();
              for (decltype(Dims) i = 0; i < Dims; ++i) {
                // As [0,0,...,0] is always contained, this operation is the same
                // as projecting onto the dimension i.
                auto PotentiallySmaller = DeltaSet.fix_si(isl::dim::set, i, 0);
                if (!PotentiallySmaller.is_equal(DeltaSet)) {
                  auto MultipleAccessSchedule = BaseSchedule.equate(
                      isl::dim::in, i, isl::dim::out, 0);
                  if (!MultipleAccessSchedule.is_equal(AlternatedLoop)) {
                    MultipleAccessSchedules.push_back(MultipleAccessSchedule);
                  }
                }
              }

              if (!MultipleAccessSchedules.empty()) {
                isl::schedule_node InnermostAccess;
                isl::schedule_node BottomUpNode = ManagedNode;
                bool FoundInnermostAccess = false;
                while (!BottomUpNode.is_equal(FilterNode)) {
                  if (isNodeType(BottomUpNode, isl_schedule_node_band)) {
                    auto Sched = bandGetPartialSchedule(BottomUpNode);
                    auto Map = filteredMultiUnionPwAffToMap(Sched, StmtId);
                    if (involveSameDimensions(AlternatedLoop, Map.remove_divs(), isl::dim::in)) {
                      FoundInnermostAccess = true;
                      InnermostAccess = BottomUpNode;
                      break;
                    }
                  }
                  BottomUpNode = BottomUpNode.parent();
                }
                auto Coincident = nodeGetCoincident(InnermostAccess, 0);
                auto CoincidentNewVal = isCoincidentLoop(Validity, AlternatedLoop);
                if (CoincidentNewVal == Coincident) {
                  LLVM_DEBUG(dbgs() << "CoincidentCheck passed " << Coincident << "\n");
                } else {
                  LLVM_DEBUG(dbgs() << "CoincidentCheck failed " << Coincident
                                    << " vs. " << CoincidentNewVal << "\n");
                  LLVM_DEBUG({
                               dbgs() << "Validity: " << Validity.to_str() << "\n";
                               dbgs() << "AlternatedLoop: " << AlternatedLoop.to_str() << "\n";
                             });
                }
                // TODO: compare coincident here;
                if (FoundInnermostAccess && Coincident) {
                  LLVM_DEBUG(dbgs() << "Innermost access:\n" << InnermostAccess.to_str() << "\n");
                  BottomUpNode = InnermostAccess.parent();
                  isl::schedule_node AlternationBaseNode;
                  isl::map AlternationBaseMap;
                  bool FoundAlternationBase = false;
                  while (!isNodeType(BottomUpNode, isl_schedule_node_domain)
                      && !FoundAlternationBase) {
                    if (isNodeType(BottomUpNode, isl_schedule_node_band)) {
                      auto Sched = bandGetPartialSchedule(BottomUpNode);
                      auto Map = filteredMultiUnionPwAffToMap(Sched, StmtId);
                      //if (Map.is_subset(AlternationBaseLoop)) {
                      //  break;
                      //}
                      // TODO: change condition
                      for (auto &PossibleSched: MultipleAccessSchedules) {
                        if (involveSameDimensions(PossibleSched, Map.remove_divs(), isl::dim::in)) {
                          FoundAlternationBase = true;
                          AlternationBaseNode = BottomUpNode;
                          AlternationBaseMap = PossibleSched;
                        }
                      }
                    }
                    BottomUpNode = BottomUpNode.parent();
                  }
                  auto AstOpts = nodeGetAstOpts(InnermostAccess);
                  auto AlternateLoopOpt = createAlternationBaseAstOption(
                      AlternationBaseMap);
                  auto PossibleOpt = AstOpts.extract_set(
                      AlternateLoopOpt.get_space());

                  if (FoundAlternationBase
                  && !PossibleOpt.is_equal(AlternateLoopOpt)) {
                    LLVM_DEBUG(dbgs() << "Shifts optimization\n");
                    auto AlternatingSched = filteredMultiUnionPwAffToMap(
                        bandGetPartialSchedule(AlternationBaseNode), StmtId);
                    auto InnermostSched = filteredMultiUnionPwAffToMap(
                        bandGetPartialSchedule(InnermostAccess), StmtId);
                    auto ScheduleEven = createAlternatingPartialSchedule(
                        AlternatingSched, InnermostSched, 0, 1);
                    auto ScheduleOdd = createAlternatingPartialSchedule(
                        AlternatingSched, InnermostSched, 1, -1);
                    auto NewPartialSchedule = ScheduleEven.unite(ScheduleOdd);
                    auto NewPartialScheduleAff
                        = isl::multi_union_pw_aff::from_union_map(NewPartialSchedule);

                    auto Permutable = nodeGetPermutable(InnermostAccess);
                    auto Coincident = nodeGetCoincident(InnermostAccess, 0);
                    AstOpts = AstOpts.unite(AlternateLoopOpt);

                    InnermostAccess
                        = isl::manage(
                        isl_schedule_node_delete(InnermostAccess.release()));
                    InnermostAccess
                        = InnermostAccess.insert_partial_schedule(NewPartialScheduleAff);
                    InnermostAccess
                        = InnermostAccess.band_member_set_coincident(0, Coincident);
                    InnermostAccess = nodeSetPermutable(InnermostAccess, Permutable);
                    InnermostAccess
                        = InnermostAccess.band_set_ast_build_options(AstOpts);

                    // search for leaf and replace
                    while (!isNodeType(InnermostAccess, isl_schedule_node_leaf)) {
                      InnermostAccess = InnermostAccess.get_child(0);
                    }
                    ManagedNode = InnermostAccess;
                    FilterNode = ManagedNode;
                    while (!isNodeType(FilterNode, isl_schedule_node_filter)
                        && !isNodeType(FilterNode, isl_schedule_node_domain)) {
                      FilterNode = FilterNode.parent();
                    }
                  }
                }
              }

              LLVM_DEBUG({
                           dbgs() << "OnlyInnerLoopAccessMap: " << OnlyInnerLoopAccessMap.to_str()
                                  << "\nAccessLexMin: " << AccessLexMin.to_str()
                                  << "\nBaseScheduleSpace: " << BaseScheduleSpace.to_str()
                                  << "\nDeltaSet: " << DeltaSet.to_str()
                                  << "\nCollected: " << Collected.to_str()
                                  << "\nApplied: " << Reversed.apply_domain(AccessLexMin).to_str()
                                  << "\nAlternated loop: " << AlternatedLoop.to_str();
                           dbgs() << "\nMultipleAccessSchedules: ";
                           for (auto &Sched: MultipleAccessSchedules) {
                             dbgs() << "\n  " << Sched.to_str();
                           }
                           dbgs() << "\n\n";
                         });
            }
          }
        }
      }
    }


    Node = ManagedNode.release();
  }
  return Node;
}

static void unifiedPatternSearch(Scop &S, const Dependences &D) {
  isl::schedule Schedule = S.getScheduleTree();
  Schedule = walkSchedule(D, Schedule, optimizeUnifiedPattern);
  LLVM_DEBUG(dbgs() << "\n\nPattern search schedule:\n" << Schedule.get_root().to_str());
  S.setScheduleTree(Schedule);
  S.markAsOptimized();
}

static void runPluto(Scop &S, const Dependences &D) {
  int ValidityKinds = Dependences::TYPE_RAW | Dependences::TYPE_WAR | Dependences::TYPE_WAW;
  isl::union_set Domain = S.getDomains();
  isl::union_map Validity = D.getDependences(ValidityKinds);
  // TODO: I have to check if this makes any difference. But I do not think so.
  Validity = Validity.gist_domain(Domain);
  Validity = Validity.gist_range(Domain);
  LLVM_DEBUG(dbgs() << "\nPluto deps: " << Validity.to_str() << "\n\n");

  // Simplify dependences - commented out to try whether Pluto can deal with
  // alternating schedules. It can not.
  //Validity = Validity.gist_domain(Domain);
  //Validity = Validity.gist_range(Domain);

  isl_ctx *Ctx = S.getIslCtx().get();

//  // set outer coincident to one to make the inner most loops more likely to be
//  // interchangeable
//  isl_options_set_schedule_outer_coincidence(Ctx, 1);
//  isl_options_set_schedule_serialize_sccs(Ctx, 1);
//  isl_options_set_schedule_maximize_band_depth(Ctx, 1);
  // 2 | 32 | 64 | 128 | 256 | 512
  auto POpts = HMPollyPlutoOptions.getValue();
  if(POpts & 1) {
    isl_options_set_schedule_serialize_sccs(Ctx, 1);
  } else {
    isl_options_set_schedule_serialize_sccs(Ctx, 0);
  }
  if(POpts & 2) {
    isl_options_set_schedule_whole_component(Ctx, 1);
  } else {
    isl_options_set_schedule_whole_component(Ctx, 0);
  }
  if(POpts & 4) {
    isl_options_set_schedule_maximize_band_depth(Ctx, 1);
  } else {
    isl_options_set_schedule_maximize_band_depth(Ctx, 0);
  }
  if(POpts & 8) {
    isl_options_set_schedule_maximize_coincidence(Ctx, 1);
  } else {
    isl_options_set_schedule_maximize_coincidence(Ctx, 0);
  }
  if(POpts & 16) {
    isl_options_set_schedule_outer_coincidence(Ctx, 1);
  } else {
    isl_options_set_schedule_outer_coincidence(Ctx, 0);
  }
  if(POpts & 32) {
    isl_options_set_schedule_algorithm(Ctx, ISL_SCHEDULE_ALGORITHM_ISL);
  } else {
    isl_options_set_schedule_algorithm(Ctx, ISL_SCHEDULE_ALGORITHM_FEAUTRIER);
  }
  if(POpts & 64) {
    isl_options_set_schedule_split_scaled(Ctx, 1);
  } else {
    isl_options_set_schedule_split_scaled(Ctx, 0);
  }
  if(POpts & 128) {
    isl_options_set_schedule_treat_coalescing(Ctx, 1);
  } else {
    isl_options_set_schedule_treat_coalescing(Ctx, 0);
  }
  if(POpts & 256) {
    isl_options_set_schedule_carry_self_first(Ctx, 1);
  } else {
    isl_options_set_schedule_carry_self_first(Ctx, 0);
  }
  if(POpts & 512) {
    isl_options_set_schedule_separate_components(Ctx, 1);
  } else {
    isl_options_set_schedule_separate_components(Ctx, 0);
  }
  isl_options_set_schedule_max_constant_term(Ctx, 20);
  isl_options_set_schedule_max_coefficient(Ctx, 20);
  isl_options_set_tile_scale_tile_loops(Ctx, 0);

  auto SC = isl::schedule_constraints::on_domain(Domain);
  SC = SC.set_validity(Validity);
  // This ensures that loops that can be alternated are marked as coincident, because
  // then they are coincident with all the dependences.
  SC = SC.set_coincidence(Validity);
  // first check with simplified deps, then do with all
  auto OldErrorSetting = isl_options_get_on_error(Ctx);
  isl_options_set_on_error(Ctx, ISL_ON_ERROR_CONTINUE);
  isl_ctx_reset_error(Ctx);
  auto Schedule = SC.compute_schedule();
  isl_error LastError = isl_ctx_last_error(Ctx);
  if(LastError) {
    isl_ctx_reset_error(Ctx);
    isl_options_set_on_error(Ctx, ISL_ON_ERROR_ABORT);
    isl::union_map Validity = D.getDependences(ValidityKinds);
    SC = isl::schedule_constraints::on_domain(Domain);
    SC = SC.set_validity(Validity);
    SC = SC.set_coincidence(Validity);
    Schedule = SC.compute_schedule();
  }
  isl_options_set_on_error(Ctx, OldErrorSetting);
  S.setScheduleTree(Schedule);

}


static void expandBandNodes(Scop &S, const Dependences &D) {
  auto Schedule = S.getScheduleTree();
  Schedule = walkSchedule(D, Schedule, splitBandNodes);
  S.setScheduleTree(Schedule);
}

static void markAsCoincident(Scop &S, const Dependences &D) {
  auto Schedule = S.getScheduleTree();
  Schedule = walkSchedule(D, Schedule, markAsCoincident);
  S.setScheduleTree(Schedule);
}

char HMRacetrackSchedulePatternMatcher::ID = 0;


static std::string getLLVMRevision() {
#ifdef LLVM_REVISION
  return LLVM_REVISION;
#else
  return "";
#endif
}

bool HMRacetrackSchedulePatternMatcher::runOnScop(Scop &S) {
  LLVM_DEBUG(dbgs() << "In racetrack schedule pattern matcher!\n");
  LLVM_DEBUG(dbgs() << "LLVM version: " << getLLVMRevision() << "\n");
  LLVM_DEBUG(dbgs() << "Shifts optimization for " << S.getNameStr() << "\n");
  LLVM_DEBUG(S.print(dbgs(), true));

  polly::hmrtm::ShiftConfig ShiftConfig(64);
  ShiftConfig.outputDirectory(HMPollyTraceDir)
  .kernelName(HMPollyKernelName)
  .traceFile(!HMPollyKernelName.empty() && !HMPollyTraceDir.empty());

  LLVM_DEBUG(ShiftConfig.print(dbgs()));
  isl_options_set_on_error(S.getIslCtx().get(), ISL_ON_ERROR_ABORT);

  LLVM_DEBUG({

    auto operations = isl_ctx_get_max_operations(S.getIslCtx().get());
    dbgs() << "Max operations: " << operations << "\n";
  });


//  auto shifts = polly::hmrtm::calculateNumberOfShifts(S, S.getSchedule(),
//      ShiftConfig.position("before-optimization"));
//  //auto shifts = 0;
//  LLVM_DEBUG(dbgs() << "\nShifts initially: " << shifts << "\n");

  LLVM_DEBUG(dbgs() << "\n\nOriginal schedule:\n" << S.getScheduleTree().get_root().to_str() << "\n\n");

  LLVM_DEBUG({
    printScheduleToStream(S.getScheduleTree(), dbgs());
  });
  const Dependences &D = getAnalysis<DependenceInfo>().getDependences(Dependences::AL_Statement);

  printAstForScop(S);


  if(D.getSharedIslCtx() != S.getSharedIslCtx()) {
    LLVM_DEBUG(dbgs() << "DependenceInfo for another SCoP/isl_ctx\n");
    return false;
  }
  if(!D.hasValidDependences()) {
    auto* Ctx = S.getIslCtx().get();
    LLVM_DEBUG(dbgs() << "No valid dependences! Recomputing...\n");
    const Dependences &D2 = getAnalysis<DependenceInfo>().recomputeDependences(Dependences::AL_Statement);
    LLVM_DEBUG(dbgs() << "successfully recomputed!\n");
    //LLVM_DEBUG(dbgs() << "f:" << isl_ctx_last_error_file(Ctx) << " l:" << isl_ctx_last_error_line(Ctx) << " m:" << isl_ctx_last_error_msg(Ctx) << "\n");
    int ValidityKinds = Dependences::TYPE_RAW | Dependences::TYPE_WAR | Dependences::TYPE_WAW;
    auto Deps = D2.getDependences(ValidityKinds);
    LLVM_DEBUG(dbgs() << "Deps: " << Deps.to_str() << "\n" << D2.hasValidDependences() << "\n");
    return false;
  }
  if(HMPollyEnablePlutoBeforePatterns) {
    LLVM_DEBUG(dbgs() << "Running Pluto...\n");
    runPluto(S, D);
    // This already calculates all permutations.
//    if(HMPollyEnableFullPermutationSearch) {
//      runAllPermutations(S, D);
//    }
    LLVM_DEBUG(dbgs() << "Expanding band nodes...\n");
//    shifts = polly::hmrtm::calculateNumberOfShifts(S, S.getSchedule(),
//        polly::hmrtm::ShiftConfig(64));
//    LLVM_DEBUG(dbgs() << "\nShifts a. Pluto:  " << shifts << "\n");
    //LLVM_DEBUG(dbgs() << "\n\nPluto schedule:\n" << S.getScheduleTree().get_root().to_str() << "\n\n");
    expandBandNodes(S, D);
    //visitTree(S.getScheduleTree());
  } else {
    LLVM_DEBUG(dbgs() << "Calculating coincident bands...\n");
    expandBandNodes(S, D);
    markAsCoincident(S, D);
  }

  LLVM_DEBUG(dbgs() << "Searching for patterns...\n");
  //LLVM_DEBUG(dbgs() << "\n\nExpanded schedule:\n" << S.getScheduleTree().get_root().to_str() << "\n\n");
  if(HMPollyEnablePatternMatching) {
    unifiedPatternSearch(S, D);
  }
//  shifts = polly::hmrtm::calculateNumberOfShifts(S, S.getSchedule(),
//      ShiftConfig.position("after-optimization"));
//  LLVM_DEBUG(dbgs() << "\nShifts a. Alter.: " << shifts << "\n");
  LLVM_DEBUG(dbgs() << "\nSchedule map: " << S.getSchedule().to_str() << "\n");
  LLVM_DEBUG(dbgs() << "\nSchedule: " << S.getScheduleTree().to_str() << "\n");

  //LLVM_DEBUG(dbgs() << "\nSchedule as UnionMap:\n" << S.getSchedule().to_str() << "\n");
  LLVM_DEBUG(dbgs() << "\n\nGenerated code:\n");

//  // Insert copy statement here and see what happens
//  // I need to create the array, a new copy statement, and the extension node
//  // with the schedule
//  // Extract the first memory access:
//  auto Stmt = S.begin();
//  auto Acc = *Stmt->begin();
//
//
//  // First the array:
//  auto ArrInfo = S.createScopArrayInfo(Acc->getElementType(), "TestCopy", {4,4});
//  // Then the statement. The tuple id of the domain and the accesses is adjusted,
//  // however, they still need to match as the domain of the memory accesses is validated.
//  auto ExtMap = isl::map(S.getIslCtx(), "{[] -> CopyStmt[i, j]: 0 <= i <= 3 and 0 <= j <= 3}");
//  auto Domain = ExtMap.range();
//  Domain = Domain.set_tuple_id(Acc->getAccessRelation().get_tuple_id(isl::dim::in));
//  auto NewAccessRel = Acc->getAccessRelation().set_tuple_id(isl::dim::out, ArrInfo->getBasePtrId());
//  auto *NewStmt = S.addScopStmt(Acc->getAccessRelation(), NewAccessRel, Domain);
//  auto ScheduleTree = S.getScheduleTree();
//  auto Node = ScheduleTree.get_root().first_child();
//  ExtMap = ExtMap.set_tuple_id(isl::dim::out, NewStmt->getDomainId());
//  auto Extension = isl::union_map(ExtMap);
//  auto NewNode = isl::schedule_node::from_extension(Extension);
//  auto ExtSched = isl::map(S.getIslCtx(), "{CopyStmt[i, j] -> [i, j]}");
//  ExtSched = ExtSched.set_tuple_id(isl::dim::in, NewStmt->getDomainId());
//
//  //NewNode = NewNode.insert_partial_schedule(isl::multi_union_pw_aff::from_union_map(isl::union_map(ExtSched))).parent();
//  LLVM_DEBUG(dbgs() << "\n" << NewNode.to_str() << "\n");
//  NewNode = NewNode.child(0).insert_partial_schedule(isl::multi_union_pw_aff::from_union_map(isl::union_map(ExtSched))).parent();
//  LLVM_DEBUG(dbgs() << "\n" << NewNode.to_str() << "\n");
//  Node = Node.graft_before(NewNode);
//  // Insert_partial_schedule only works if you use the (hidden) leaf which is
//  // a child of the extension node as it inserts the schedule above the node.
//  LLVM_DEBUG({
//    dbgs() << "\n\n" << Node.to_str() << "\n\n";
//  });
//  ScheduleTree = Node.get_schedule();
//
//  ScheduleTree = hoistExtensionNodes(Node.get_schedule());
//  LLVM_DEBUG(dbgs() << "\n\nNew schedule: " << ScheduleTree.get_root().to_str());
//  S.setScheduleTree(ScheduleTree);
//  LLVM_DEBUG({
//               dbgs() << "\n\n\n";
//               S.print(dbgs(), true);
//             });
  printAstForScop(S);

  //const Dependences &D2 = getAnalysis<DependenceInfo>().recomputeDependences(Dependences::AL_Statement);

  //runPluto(S, D2);
//  shifts = polly::hmrtm::calculateNumberOfShifts(S, S.getSchedule());
//  LLVM_DEBUG(dbgs() << "\nShifts after:  " << shifts << "\n");

  // This code is based on the normal polly optimizer pass, except that
  // it performs only the alternation optimization after the Pluto run.
  return false;
}

void HMRacetrackSchedulePatternMatcher::getAnalysisUsage(AnalysisUsage &AU) const {
  ScopPass::getAnalysisUsage(AU);
  AU.addRequired<DependenceInfo>();
  // Dependence info is NOT preserved by this pass...
  //AU.addPreserved<DependenceInfo>();
}

Pass *polly::createHMRacetrackSchedulePatternMatcherPass() {
  return new HMRacetrackSchedulePatternMatcher();
}

/// Unused Permutation code for now.

static isl_stat moveToNextLexicographic(isl::schedule_node &Node) {
  // always go to the left most child first
  if(Node.has_children()) {
    Node = Node.first_child();
    return isl_stat_ok;
  }
  while(Node.has_parent() && !Node.has_next_sibling()) {
    Node = Node.parent();
  }
  if(Node.has_next_sibling()) {
    Node = Node.next_sibling();
    return isl_stat_ok;
  }
  return isl_stat_error;
}

static isl_stat moveToNextPermutableBand(isl::schedule_node &Node) {
  while(!isNodeType(Node, isl_schedule_node_band) || !nodeGetPermutable(Node)
      || isl_schedule_node_band_n_member(Node.get()) < 2) {
    auto stat = moveToNextLexicographic(Node);
    if(stat == isl_stat_error) {
      return isl_stat_error;
    }
  }
  return isl_stat_ok;
}

template<typename Callback>
static void generatePermutations(isl::schedule_node &Node, Callback &Cb);

template<typename Callback>
static void outputPermutation(isl::schedule_node Node, Callback &Cb) {
  auto Stat = moveToNextLexicographic(Node);
  if(Stat != isl_stat_error) {
    Stat = moveToNextPermutableBand(Node);
  }
  // finished permuting, reached end of tree
  if(Stat == isl_stat_error) {
    LLVM_DEBUG(dbgs() << "Outputting permutation!\n");
    Cb(Node.get_schedule());
  } else {
    // create permutations for next band
    generatePermutations(Node, Cb);
  }
}

template<typename Callback>
static void generatePermutations(isl::schedule_node &Node, Callback &Cb) {
  // implementing the heap algorith here, iterative version
  //https://en.wikipedia.org/wiki/Heap%27s_algorithm
  auto BandMembers = nodeBandMembers(Node);
  std::vector<int> StackState(BandMembers); // c from wiki
  std::fill(StackState.begin(), StackState.end(), 0);
  outputPermutation(Node, Cb);
  // stack pointer
  int i = 0;
  while(i < BandMembers) {
    if(StackState[i] < i) {
      if(i % 2 == 0) {
        Node = permuteBandNodeDimensions(Node, 0, i);
      } else {
        Node = permuteBandNodeDimensions(Node, StackState[i], i);
      }
      outputPermutation(Node, Cb);
      ++StackState[i];
      i = 0;
    } else {
      StackState[i] = 0;
      ++i;
    }
  }
}

template<typename Callback>
static void visitTree(const isl::schedule &Tree, Callback cb ) {
  isl::schedule_node Node = Tree.get_root();
  do {
    LLVM_DEBUG(dbgs() << Node.to_str() << "\n\n");
  }while(moveToNextLexicographic(Node) != isl_stat_error);
}
static void runAllPermutations(Scop &S, const Dependences &D) {
  LLVM_DEBUG(dbgs() << "\n\nPERMUTATION RUNNER!\n");
  auto Lambda = [&D, &S](isl::schedule Schedule) {
    Schedule = walkSchedule(D, Schedule, splitBandNodes);
    Schedule = walkSchedule(D, Schedule, optimizeUnifiedPattern);
    auto shifts = polly::hmrtm::calculateNumberOfShifts(S, Schedule.get_map());
    LLVM_DEBUG(dbgs() << "Shifts in enumeration:" << shifts << "\n"
                      << Schedule.get_root().to_str() << "\n");
  };
  outputPermutation(S.getScheduleTree().get_root(), Lambda);
  LLVM_DEBUG(dbgs() << "END PERMUTATION RUNNER!\n\n");
}




INITIALIZE_PASS_BEGIN(HMRacetrackSchedulePatternMatcher,
    "polly-opt-hm-pattern",
    "Polly Racetrack Optimizer - Custom pass by Hauke Mewes", false, false)
INITIALIZE_PASS_DEPENDENCY(DependenceInfo)
INITIALIZE_PASS_DEPENDENCY(ScopInfoRegionPass)
INITIALIZE_PASS_END(HMRacetrackSchedulePatternMatcher,
    "polly-opt-hm-pattern",
    "Polly Racetrack Optimizer - Custom pass by Hauke Mewes", false, false)

#undef DEBUG_TYPE