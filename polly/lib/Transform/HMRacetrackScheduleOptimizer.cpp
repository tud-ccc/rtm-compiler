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
#include "polly/DependenceInfo.h"
#include "polly/LinkAllPasses.h"
#include "polly/ScheduleOptimizer.h"
#include "polly/ScopInfo.h"

#include "llvm/Support/CommandLine.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/raw_ostream.h"

#include "isl/ctx.h"
#include "isl/options.h"
#include "isl/printer.h"
#include "isl/schedule.h"
#include "isl/schedule_node.h"
#include "isl/union_map.h"
#include "isl/union_set.h"

#define DEBUG_TYPE "polly-opt-hm-racetrack"


using namespace llvm;
using namespace polly;

namespace {
class HMRacetrackScheduleOptimizer : public ScopPass {
public:
  static char ID;

  explicit HMRacetrackScheduleOptimizer() : ScopPass(ID) {}

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
isl::schedule walkSchedule(const Dependences &D, isl::schedule &Schedule,
    __isl_give isl_schedule_node *(*fn)(
        __isl_take isl_schedule_node *, void *)) {
  auto Root = Schedule.get_root();
  Root = isl::manage(isl_schedule_node_map_descendant_bottom_up(
      Root.release(), fn, const_cast<void *>(static_cast<const void *>(&D))
  ));
  return Root.get_schedule();
}

/**
 * This creates the isl map that can alternate the inner most loop depending
 * on the sign of the former two.
 * For example, with space size 3, modulus constant 2 and output sign -1, it
 * creates:
 * { Stmt3[i0, i1, i2] -> [o0] :
 * exists (e0 = floor((1 + i0 + i1)/2): o0 = -i2 and 2e0 = 1 + i0 + i1) }
 * @param Ctx Isl context.
 * @param spaceSize Total number of dimensions.
 * @param spaceIdentifier Statement identifier.
 * @param modulusConstant Modulus on which to apply this. Usually equals 2.
 * @param outputSign Should the iteration go forward or backward?
 * @return The appropriate partial map.
 */
static isl::map createAffineForModulus(isl::ctx Ctx, unsigned spaceSize, isl::id spaceIdentifier, int modulusConstant, int outputSign) {
  // This is a good example of isl. To create a map, first the required space is
  // created.
  isl::space space(Ctx, 0, spaceSize, spaceSize);
  // Then the identifier of the space has to be set.
  space = space.set_tuple_id(isl::dim::in, spaceIdentifier);
  // To create a map that one wants, it is best to start with the universe map.
  isl::map affAsMap = isl::map::universe(space);
  // Drop all dimensions but two from the output space.
  affAsMap = affAsMap.project_out(isl::dim::out, 0, spaceSize - 2);
  // To create a constraint, one has to create a local space from the space.
  isl::local_space localSpace(affAsMap.get_space());
  // Create modulus constraint on second last input input.
  isl::constraint modulusConstraint = isl::constraint::alloc_equality(localSpace);
  // Adjust coefficients
  modulusConstraint = modulusConstraint.set_coefficient_si(isl::dim::in, spaceSize - 2, 1);
  modulusConstraint = modulusConstraint.set_coefficient_si(isl::dim::in, spaceSize - 3, 1);
  modulusConstraint = modulusConstraint.set_coefficient_si(isl::dim::out, 0, -2);
  modulusConstraint = modulusConstraint.set_constant_si(modulusConstant);
  affAsMap = affAsMap.add_constraint(modulusConstraint);
  // Add the output constraint.
  isl::constraint outputConstraint = isl::constraint::alloc_equality(localSpace);
  outputConstraint = outputConstraint.set_coefficient_si(isl::dim::in, spaceSize - 1, outputSign);
  outputConstraint = outputConstraint.set_coefficient_si(isl::dim::out, 1, -1);
  affAsMap = affAsMap.add_constraint(outputConstraint);
  // Make the modulus constraint an existential quantifier.
  affAsMap = affAsMap.project_out(isl::dim::out, 0, 1);
  return affAsMap;
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
/// Restore the initial ordering of dimensions of the band node
///
/// In case the band node represents all the dimensions of the iteration
/// domain, recreate the band node to restore the initial ordering of the
/// dimensions.
///
/// @param Node The band node to be modified.
/// @return The modified schedule node.
static isl::schedule_node
getBandNodeWithOriginDimOrder(isl::schedule_node Node) {
  assert(isl_schedule_node_get_type(Node.get()) == isl_schedule_node_band);
  if (isl_schedule_node_get_type(Node.child(0).get()) != isl_schedule_node_leaf)
    return Node;
  auto Domain = Node.get_universe_domain();
  assert(isl_union_set_n_set(Domain.get()) == 1);
  if (Node.get_schedule_depth() != 0 ||
      (isl::set(Domain).dim(isl::dim::set) !=
          isl_schedule_node_band_n_member(Node.get())))
    return Node;
  Node = isl::manage(isl_schedule_node_delete(Node.copy()));
  auto PartialSchedulePwAff = Domain.identity_union_pw_multi_aff();
  auto PartialScheduleMultiPwAff =
      isl::multi_union_pw_aff(PartialSchedulePwAff);
  PartialScheduleMultiPwAff =
      PartialScheduleMultiPwAff.reset_tuple_id(isl::dim::set);
  return Node.insert_partial_schedule(PartialScheduleMultiPwAff);
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

/**
 * Applies the alternating loop transformation to the schedule if the MMA is
 * detected.
 * @param S The scop to analyze and modify.
 * @param D The dependence analysis of the scop.
 * @param Schedule The current schedule for the scop.
 */
static void optimizeMatMulPattern(Scop &S, const Dependences &D,
    isl::schedule &Schedule) {
  auto optimizeBandNode = [] (__isl_take isl_schedule_node * Node, void
  * DAsVoid) -> __isl_give isl_schedule_node * {
    if (!ScheduleTreeOptimizer::isTileableBandNode(isl::manage_copy(Node))) {
      return Node;
    }
    auto D = static_cast<const Dependences *>(DAsVoid);

    MatMulInfoTy MMI;
    if (ScheduleTreeOptimizer::isMatrMultPattern(isl::manage_copy(Node), D,
        MMI)) {
      // That means this node actually contains the raw matmul pattern
      // Now I only have to apply my loop transformation again.

      LLVM_DEBUG({
        dbgs() << "Found MMA!\n";
        dbgs() << "MMI.i =" << MMI.i << "\nMMI.j = " << MMI.j << "\nMMI.k = "
          << MMI.k << "\n";
      });
      auto ManagedNode = isl::manage(Node);
      // Still taken from the original matmul optimization code.
      // This moves the three iterations to the three innermost ones.
      int DimOutNum = isl_schedule_node_band_n_member(ManagedNode.get());
      ManagedNode = permuteBandNodeDimensions(ManagedNode, MMI.i, DimOutNum - 3);
      int NewJ = MMI.j == DimOutNum - 3 ? MMI.i : MMI.j;
      int NewK = MMI.k == DimOutNum - 3 ? MMI.i : MMI.k;
      ManagedNode = permuteBandNodeDimensions(ManagedNode, NewJ, DimOutNum - 2);
      NewK = NewK == DimOutNum - 2 ? NewJ : NewK;
      ManagedNode = permuteBandNodeDimensions(ManagedNode, NewK, DimOutNum - 1);
      Node = ManagedNode.release();
      // Split the last band
      // This code does not yet work if there exist more inner loops.
      Node = isl_schedule_node_band_split(Node,
          isl_schedule_node_band_n_member(Node) - 1);
      ManagedNode = isl::manage(Node);
      ManagedNode = ManagedNode.first_child();
      // extract the first domain, should be the only one.
      isl::set Domain = ManagedNode.get_domain().get_set_list().get_at(0);
      isl::ctx Ctx = ManagedNode.get_ctx();
      auto domainDim = Domain.n_dim();
      // Create the two schedule for even and odd iterations.
      isl::union_map evenUnionAff(createAffineForModulus(Ctx, domainDim, Domain.get_tuple_id(), 0, 1));
      isl::union_map oddUnionAff(createAffineForModulus(Ctx, domainDim, Domain.get_tuple_id(), 1, -1));
      // ... and combine them.
      isl::union_map combined = evenUnionAff.unite(oddUnionAff);
      // Insert into the right place of the schedule tree.
      auto aff = isl::multi_union_pw_aff::from_union_map(combined);
      ManagedNode = ManagedNode.insert_partial_schedule(aff);
      ManagedNode = ManagedNode.get_child(0);
      ManagedNode = ManagedNode.cut();
      ManagedNode = ManagedNode.parent();
      Node = ManagedNode.release();
    }
    return Node;
  };
  // Search for the optimization.
  isl::schedule newSchedule = walkSchedule(D, Schedule, optimizeBandNode);
  S.setScheduleTree(newSchedule);
  newSchedule.dump();
}

char HMRacetrackScheduleOptimizer::ID = 0;

static __isl_give isl_schedule_node * splitBandNodes(__isl_take isl_schedule_node *Node, void *User) {
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

static isl::multi_union_pw_aff bandGetPartialSchedule(
    const isl::schedule_node &Node) {
  // TODO: what about get_partial_schedule_union_map?
  return isl::manage(isl_schedule_node_band_get_partial_schedule(Node.get()));
}

static isl_bool nodeGetCoincident(const isl::schedule_node &Node, int Pos) {
  return isl_schedule_node_band_member_get_coincident(Node.get(), Pos);
}

static isl::boolean isCoincidentLoop(const isl::union_map &Deps,
                                     const isl::union_map &Schedule) {
  auto DeltaSet = Deps.apply_domain(Schedule).apply_range(Schedule)
      .deltas().detect_equalities().coalesce();
  return DeltaSet.is_empty() || DeltaSet.is_equal(isl::set(Deps.get_ctx(), "{[0]}"));

}

static __isl_give isl_schedule_node * checkCoincident(__isl_take isl_schedule_node *Node, void *User) {
  if(isl_schedule_node_get_type(Node) == isl_schedule_node_band) {
    auto D = static_cast<Dependences *>(User);
    int ValidityKinds = Dependences::TYPE_RAW | Dependences::TYPE_WAR | Dependences::TYPE_WAW;
    isl::union_map Validity = D->getDependences(ValidityKinds);
    auto ManagedNode = isl::manage(Node);
    auto AffSched = bandGetPartialSchedule(ManagedNode);
    auto Sched = isl::union_map::from(AffSched);
    //auto SchedList = Scheds.get_map_list();
    auto Coincident = isl::manage(nodeGetCoincident(ManagedNode, 0));
    //for(int i = 0; i < SchedList.n_map(); ++i) {
      //auto Sched = SchedList.get_at(i);
      auto CalcCoincident = isCoincidentLoop(Validity, Sched);
      if(CalcCoincident == Coincident) {
        LLVM_DEBUG(dbgs() << "CoincidentCheck passed " << Coincident << "\n");
      } else {
        LLVM_DEBUG(dbgs() << "CoincidentCheck failed " << Coincident
                          << " vs. " << CalcCoincident << "\n");
        LLVM_DEBUG({
                     dbgs() << "CoincidentCheck Validity: " << Validity.to_str() << "\n";
                     dbgs() << "CoincidentCheck AlternatedLoop: " << Sched.to_str() << "\n";
                   });
      }
    //}
    Node = ManagedNode.release();
  }
  return Node;
}

bool HMRacetrackScheduleOptimizer::runOnScop(Scop &S) {
  // This code is based on the normal polly optimizer pass, except that
  // it performs only the alternation optimization after the Pluto run.

  LLVM_DEBUG(dbgs() << "In Racetrack Schedule Optimizer!\n");
  LLVM_DEBUG(dbgs() << "CoincidentCheck on: " << S.getNameStr() << "\n");
  //LLVM_DEBUG(dbgs() << "Schedule tree before Pluto:\n" << S.getScheduleTree().to_str() << "\n");
  //S.getScheduleTree().dump();
  //auto shifts = polly::hmrtm::calculateNumberOfShifts(S, S.getSchedule());
  //printAstForScop(S);
  //LLVM_DEBUG(dbgs() << "Shifts: " << shifts << "\n");

  // Skip empty SCoPs but still allow code generation as it will delete the
  // loops present but not needed.
  if (S.getSize() == 0)  {
    S.markAsOptimized();
    return false;
  }

  const Dependences &D = getAnalysis<DependenceInfo>().getDependences
      (Dependences::AL_Statement);

  if (D.getSharedIslCtx() != S.getSharedIslCtx()) {
    LLVM_DEBUG(dbgs() << "DependenceInfo for another SCoP/isl_ctx!\n");
    return false;
  }

  if (!D.hasValidDependences()) {
    return false;
  }

  int ValidityKinds = Dependences::TYPE_RAW | Dependences::TYPE_WAR |
      Dependences::TYPE_WAW;
  int ProximityKinds = ValidityKinds;

  isl::union_set Domain = S.getDomains();

  if(!Domain) {
    return false;
  }

  isl::union_map Validity = D.getDependences(ValidityKinds);
  isl::union_map Proximity = D.getDependences(ProximityKinds);

  // Simplify the dependences by removing the constraints introduced by the
  // domains. See the ScheduleOptimizer for details.
//  Validity = Validity.gist_domain(Domain);
//  Validity = Validity.gist_range(Domain);
//  Proximity = Proximity.gist_domain(Domain);
//  Proximity = Proximity.gist_range(Domain);

  isl_ctx *Ctx = S.getIslCtx().get();

  isl_options_set_schedule_outer_coincidence(Ctx, 0);
  isl_options_set_schedule_serialize_sccs(Ctx, 1);
  isl_options_set_schedule_maximize_band_depth(Ctx, 1);
  isl_options_set_schedule_max_constant_term(Ctx, 20);
  isl_options_set_schedule_max_coefficient(Ctx, 20);
  isl_options_set_tile_scale_tile_loops(Ctx, 0);

  auto SC = isl::schedule_constraints::on_domain(Domain);
  SC = SC.set_proximity(Proximity);
  SC = SC.set_validity(Validity);
  SC = SC.set_coincidence(Validity);

  auto OnErrorStatus = isl_options_get_on_error(Ctx);
  isl_options_set_on_error(Ctx, ISL_ON_ERROR_CONTINUE);
  auto Schedule = SC.compute_schedule();
  isl_options_set_on_error(Ctx, OnErrorStatus);
  if(!Schedule) {
    LLVM_DEBUG(dbgs() << "CoincidentCheck: Pluto fail!\n");
  }

  Schedule = walkSchedule(D, Schedule, splitBandNodes);
  Schedule = walkSchedule(D, Schedule, checkCoincident);

//  LLVM_DEBUG(dbgs() << "Schedule tree after Pluto:\n");
//  Schedule.dump();
//  optimizeMatMulPattern(S, D, Schedule);
//
//  S.getSchedule().dump();
//  shifts = polly::hmrtm::calculateNumberOfShifts(S, S.getSchedule());
//  LLVM_DEBUG(dbgs() << "Shifts: " << shifts << "\n");
  return false;
}

void HMRacetrackScheduleOptimizer::getAnalysisUsage(AnalysisUsage &AU) const {
  ScopPass::getAnalysisUsage(AU);
  AU.addRequired<DependenceInfo>();
  AU.addPreserved<DependenceInfo>();
}

Pass *polly::createHMRacetrackScheduleOptimizerPass() {
  return new HMRacetrackScheduleOptimizer();
}

INITIALIZE_PASS_BEGIN(HMRacetrackScheduleOptimizer,
    "polly-opt-hm-racetrack",
    "Polly Racetrack Optimizer - Custom pass by Hauke Mewes", false, false)
INITIALIZE_PASS_DEPENDENCY(DependenceInfo)
INITIALIZE_PASS_DEPENDENCY(ScopInfoRegionPass)
INITIALIZE_PASS_END(HMRacetrackScheduleOptimizer,
    "polly-opt-hm-racetrack",
    "Polly Racetrack Optimizer - Custom pass by Hauke Mewes", false, false)

#undef DEBUG_TYPE