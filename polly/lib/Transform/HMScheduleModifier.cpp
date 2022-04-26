#include "polly/DependenceInfo.h"
#include "polly/LinkAllPasses.h"
#include "polly/Options.h"
#include "polly/ScopInfo.h"
#include "polly/HMScheduleReconstruction.h"

#include "isl/isl-noexceptions.h"

#include "llvm/Analysis/DependenceAnalysis.h"
#include "llvm/Support/CommandLine.h"

using namespace llvm;
using namespace polly;
#define DEBUG_TYPE "hm-polly"

static cl::opt<bool>
    HMScheduleModifierEnabled("polly-hm-schedule-modifier-enabled", cl::desc("Enable Hauke Mewes' schedule modifier"), cl::init(true), cl::ZeroOrMore);
static cl::opt<bool>
    HMScheduleModifierAnalyzeOnly("polly-hm-schedule-modifier-analyze-only", cl::desc("Put Hauke Mewes' schedule modifier in analyze only mode."),
        cl::init(false), cl::ZeroOrMore);
// Use anonymous namespace for llvm pass.
// This is common best practice.
namespace {

class HMScheduleModifier : public ScopPass {
public:
  static char ID;
  explicit HMScheduleModifier() : ScopPass(ID) {
    LLVM_DEBUG(dbgs() << "Create schedule modifier pass!\n");
  }

  bool runOnScop(Scop &S) override;

  void fiddle(Scop &S);

  bool reorderInnnerLoopMemAccess(Scop &S);

  /// Register all analyses and transformation required.
  void getAnalysisUsage(AnalysisUsage &AU) const override;
private:
  isl::map createAffineForModulus(isl::ctx Ctx, int spaceSize, isl::id spaceIdentifier, int modulusConstant, int outputSign);
  void analyzeModifiedScheduleAndDependencies(Scop &S);
};
}


char HMScheduleModifier::ID = 0;

isl::map HMScheduleModifier::createAffineForModulus(isl::ctx Ctx, int spaceSize, isl::id spaceIdentifier, int modulusConstant, int outputSign) {
  isl::space space(Ctx, 0, spaceSize, spaceSize);
  space = space.set_tuple_id(isl::dim::in, spaceIdentifier);
  isl::map affAsMap = isl::map::universe(space);
  // drop all dimensions but two from the outer space
  affAsMap = affAsMap.project_out(isl::dim::out, 0, spaceSize - 2);
  isl::local_space localSpace(affAsMap.get_space());
  // create modulus constraint on second last input input
  isl::constraint modulusConstraint = isl::constraint::alloc_equality(localSpace);
  modulusConstraint = modulusConstraint.set_coefficient_si(isl::dim::in, spaceSize - 2, 1);
  modulusConstraint = modulusConstraint.set_coefficient_si(isl::dim::in, spaceSize - 3, 1);
  modulusConstraint = modulusConstraint.set_coefficient_si(isl::dim::out, 0, -2);
  modulusConstraint = modulusConstraint.set_constant_si(modulusConstant);
  affAsMap = affAsMap.add_constraint(modulusConstraint);
  isl::constraint outputConstraint = isl::constraint::alloc_equality(localSpace);
  outputConstraint = outputConstraint.set_coefficient_si(isl::dim::in, spaceSize - 1, outputSign);
  outputConstraint = outputConstraint.set_coefficient_si(isl::dim::out, 1, -1);
  affAsMap = affAsMap.add_constraint(outputConstraint);
  // make the modulus constraint an existential quantifier
  affAsMap = affAsMap.project_out(isl::dim::out, 0, 1);
  return affAsMap;
}


static void printAstForScop(Scop &S) {
  auto ast = isl::ast_build::from_context(S.getContext());
  auto printer = isl_printer_to_str(ast.get_ctx().get());
  printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
  auto scheduleTree = S.getScheduleTree();
  auto astNode = ast.node_from_schedule(scheduleTree);

  printer = isl_printer_print_ast_node(printer, astNode.get());
  std::string out(isl_printer_get_str(printer));
  LLVM_DEBUG(dbgs() << out << "\n");
  isl_printer_free(printer);
}

void HMScheduleModifier::analyzeModifiedScheduleAndDependencies(Scop &S) {
  auto scheduleTree = S.getScheduleTree();
  LLVM_DEBUG(dbgs() << "In analysis\n");
  auto &dInfo = getAnalysis<polly::DependenceInfo>();
  dInfo.recomputeDependences(Dependences::AL_Statement);
  const auto &d = dInfo.getDependences(Dependences::AL_Statement);
  auto RAW = d.getDependences(Dependences::TYPE_RAW);
  RAW.dump();
  auto WAR = d.getDependences(Dependences::TYPE_WAR);
  WAR.dump();
  auto WAW = d.getDependences(Dependences::TYPE_WAW);
  WAW.dump();
  auto all = d.getDependences(Dependences::TYPE_RAW | Dependences::TYPE_WAR | Dependences::TYPE_WAW);
  all.dump();
  // next step: use these dependences and a newly constructed domain with the scheduler constraints in order to see if we get
  // the same schedule again
  //scheduleTree.dump();
}

bool HMScheduleModifier::reorderInnnerLoopMemAccess(Scop &S) {
  isl::schedule ScheduleTree = S.getScheduleTree();
  isl::union_set Domains = S.getDomains();
  // get the first domain from the domains list and use its space
  isl::set Domain = Domains.get_set_list().get_at(0);
  isl::ctx Ctx = Domains.get_ctx();
  auto domainDim = Domain.n_dim();
  if(domainDim < 3) {
    return false;
  }
  isl::union_map evenUnionAff(createAffineForModulus(Ctx, domainDim, Domain.get_tuple_id(), 0, 1));
  isl::union_map oddUnionAff(createAffineForModulus(Ctx, domainDim, Domain.get_tuple_id(), 1, -1));
  isl::union_map combined = evenUnionAff.unite(oddUnionAff);
  auto aff = isl::multi_union_pw_aff::from_union_map(combined);
  //aff.dump();
  auto node = ScheduleTree.get_root();
  // assume we have exactly one band node per loop
  // descend to the inner most loop
  for(decltype(domainDim) i = 0; i < domainDim; ++i) {
    node = node.get_child(0);
  }
  //node.get_prefix_schedule_multi_union_pw_aff().dump();
  node = node.insert_partial_schedule(aff);
  node = node.get_child(0);
  node = node.cut();
  LLVM_DEBUG(dbgs() << "Modifying schedule...\n");
  node.get_schedule().dump();
  node.get_schedule().get_map().dump();
  if(!HMScheduleModifierAnalyzeOnly) {
    S.setScheduleTree(node.get_schedule());
    S.markAsOptimized();
  }
  printAstForScop(S);


  return false;
}

bool HMScheduleModifier::runOnScop(Scop &S) {
  //S.getScheduleTree().dump();
  //printAstForScop(S);

  auto tree = S.getScheduleTree();
  isl_schedule_node_foreach_descendant_top_down(tree.get_root().get(), [](isl_schedule_node * nodeptr, void * second) -> isl_bool {
    // just take a copy of this node
    isl::schedule_node node = isl::manage_copy(nodeptr);
    switch(isl_schedule_node_get_type(node.get())) {
    case isl_schedule_node_band:
      break;
    case isl_schedule_node_filter:
      LLVM_DEBUG({
        node.filter_get_filter().dump();
      });
      break;
    default:
      break;
    }

    return isl_bool_true;
  }, nullptr);

  if(!HMScheduleModifierEnabled) {
    LLVM_DEBUG(dbgs() << "Not allowed to run!\n");
    return false;
  }

  // This would allow to access the dependence analysis from Preston Briggs.
  //auto &dependence = getAnalysis<DependenceAnalysisWrapperPass>();
  //dependence.dump();
  const Dependences &D =
      getAnalysis<polly::DependenceInfo>().getDependences(Dependences::AL_Statement);
  isl::union_map raw = D.getDependences(Dependences::TYPE_RAW);
  raw.dump();
  auto accesses = S.getAccesses();
  for(auto &statement: S) {
    for(auto access: statement) {
      LLVM_DEBUG(dbgs() << "Access relation string: " << access->getAccessRelationStr() << "\n");
      access->getAddressFunction().dump();
      access->getAccessRelation().get_tuple_id(isl::dim::out).dump();
      std::string opcode(access->getAccessInstruction()->getOpcodeName());
      LLVM_DEBUG(dbgs() << "Opcode " << opcode << "\n");
    }
  }
  LLVM_DEBUG(dbgs() << "Analyzing block...\n");
  for(auto &statement: S) {
    auto &instructions = statement.getInstructions();
    for(Instruction* i: instructions) {
      LLVM_DEBUG(dbgs() << "Opcode: " << i->getOpcodeName() << "\n");
    }
  }



  bool result = reorderInnnerLoopMemAccess(S);
  //analyzeModifiedScheduleAndDependencies(S);
  return result;

//  LLVM_DEBUG(dbgs() << "Test debug flag!\n");
//  isl::union_map Schedule = S.getSchedule();
//  isl::union_set Domain = S.getDomains();
//  isl::schedule ScheduleTree = S.getScheduleTree();
//  //ScheduleTree.dump();
//
//  // try to investigate the set list
//  auto setList = Domain.get_set_list();
//  isl::id tupleId;
//  isl::space space;
//  for(auto i = 0; i < setList.size(); ++i) {
//    auto set = setList.get_at(i);
//    set.get_space().dump();
//    space = set.get_space();
//    tupleId = set.get_tuple_id();
//    auto tupleName = set.get_tuple_name();
//    LLVM_DEBUG(dbgs() <<"In domain set:" << i << "\n");
//    LLVM_DEBUG(dbgs() << "Tuple name: " << tupleName << "\n");
//    tupleId.dump();
//  }
//
////  Schedule.dump();
////  Domain.dump();
//  // this is executed directly after the potential isl opt pass
//  auto Root = ScheduleTree.get_root();
//  auto FirstChild = Root.get_child(0);
//  auto LevelTwoChild = FirstChild.get_child(0);
//  auto LevelThreeChild = LevelTwoChild.get_child(0);
//
//  //auto Space = isl::manage(isl_schedule_node_band_get_space(LevelThreeChild.get()));
//  //auto Dims = Space.dim(isl::dim::set);
//  //auto Sizes = isl::multi_val::zero(Space);
//  /*auto isBand = isl_schedule_node_get_type(FirstChild.get()) == isl_schedule_node_band;
//  LLVM_DEBUG(dbgs() << "Is band?" << (isBand ? "yes" : "no") << "\n");
//  LLVM_DEBUG(dbgs() << "Members?" << isl_schedule_node_band_n_member(FirstChild.get()) << "\n");*/
////  auto prefixScheduleFirstChild = FirstChild.get_prefix_schedule_multi_union_pw_aff();
////  auto prefixScheduleLevelTwoChild = LevelTwoChild.get_prefix_schedule_multi_union_pw_aff();
////  auto prefixScheduleLevelThreeChild = LevelThreeChild.get_prefix_schedule_multi_union_pw_aff();
////
////  auto unionMap = FirstChild.get_subtree_schedule_union_map();
////  unionMap.dump();
////  LevelTwoChild.get_subtree_schedule_union_map().dump();
////  LevelThreeChild.get_subtree_schedule_union_map().dump();
//
//
//  LLVM_DEBUG(dbgs() << "Begin explore\n");
//  auto prefixSchedule = LevelThreeChild.get_prefix_schedule_union_pw_multi_aff();
//  auto list = prefixSchedule.get_pw_multi_aff_list();
//  for(int i = 0; i < list.size(); ++i) {
//    auto part = list.get_at(i);
//    part.dump();
//    auto space = part.get_space();
//    space.dump();
//    LLVM_DEBUG(dbgs() << space.dim(isl::dim::in) << "\n");
//    LLVM_DEBUG(dbgs() << space.dim(isl::dim::out) << "\n");
//  }
//  LLVM_DEBUG(dbgs() << "End explore\n");
//
//  // list dependences
//  const Dependences &D =
//      getAnalysis<DependenceInfo>().getDependences(Dependences::AL_Statement);
//  isl::union_map raw = D.getDependences(Dependences::TYPE_RAW);
//  raw.get_space().dump();
//  raw = raw.gist_domain(Domain);
//  isl::union_map war = D.getDependences(Dependences::TYPE_WAR);
//  war = war.gist_range(Domain);
//  isl::union_map waw = D.getDependences(Dependences::TYPE_WAW);
//  LLVM_DEBUG(dbgs() << "Domain: \n");
//  Domain.dump();
//  LLVM_DEBUG(dbgs() << "RAW: \n");
//  raw.dump();
//  LLVM_DEBUG(dbgs() << "WAR: \n");
//  war.dump();
//  LLVM_DEBUG(dbgs() << "WAW: \n");
//  waw.dump();
//
//  auto ctx = S.getIslCtx();
//  auto validityEven = isl::map(ctx, "{ Stmt2[i0, i1, i2] -> Stmt2[i0, i1, 1 + i2] : (i0 + i1) % 2 = 0}");
//  validityEven.dump();
//  auto validityOdd = isl::map(ctx, "{ Stmt2[i0, i1, i2] -> Stmt2[i0, i1, i2 - 1]: (i0 + i1) % 2 = 1}");
//
//  validityOdd.dump();
//  validityOdd = validityOdd.set_tuple_id(isl::dim::in, tupleId);
//  validityOdd = validityOdd.set_tuple_id(isl::dim::out, tupleId);
//  validityEven = validityEven.set_tuple_id(isl::dim::in, tupleId);
//  validityEven = validityEven.set_tuple_id(isl::dim::out, tupleId);
//
//  LLVM_DEBUG(dbgs() << "dimensions in: " <<validityEven.dim(isl::dim::in) << "\n");
//  validityEven.get_dim_id(isl::dim::in, 0).dump();
//  validityEven.get_space().dump();
////  auto mapList = validityOdd.get_map_list();
////  LLVM_DEBUG(dbgs() << "mapList size: " << mapList.size() << "\n");
////  for(auto i = 0; i < mapList.size(); ++i) {
////    auto map = mapList.get_at(i);
////    auto tupleNameIn = map.get_tuple_name(isl::dim::in);
////    auto tupleIdIn = map.get_tuple_id(isl::dim::in);
////    tupleIdIn.dump();
////  }
//
//  isl::local_space localSpace(space);
//  isl::constraint local = isl::constraint::alloc_inequality(localSpace);
//  local.dump();
//
//  auto unionValidities = validityEven.unite(validityOdd);
//  unionValidities.dump();
//
//  auto sc = isl::schedule_constraints::on_domain(Domain);
//  sc = sc.set_validity(validityEven);
//  //sc = sc.set_coincidence(validityOdd);
//  auto schedule2 = sc.compute_schedule();
//  schedule2.dump();
//
//
//  // only transform if analyzer mode is not active.
//  if(!HMScheduleModifierAnalyzeOnly) {
//    LLVM_DEBUG(dbgs() << "Modifying schedule...\n");
//    FirstChild = isl::manage(isl_schedule_node_band_sink(FirstChild.release()));
//    LevelTwoChild = FirstChild.get_child(0);
//    LevelTwoChild = isl::manage(isl_schedule_node_band_sink(LevelTwoChild.release()));
//    S.setScheduleTree(LevelTwoChild.get_schedule());
//    S.markAsOptimized();
//  }
//
//
//  fiddle(S);
//
//  //prefixScheduleFirstChild.dump();
//  //prefixScheduleLevelTwoChild.dump();
//  //prefixScheduleLevelThreeChild.dump();
//
//  //S.getScheduleTree().dump();
//
//  return false;
}

void HMScheduleModifier::fiddle(Scop &S) {
  LLVM_DEBUG(dbgs() << "In fiddle...\n");
  auto ScopDomains = S.getDomains();
  auto ScopSchedule = S.getSchedule();
  auto ScopDomainSet = ScopDomains.get_set_list().get_at(0);
  auto ScopSpace = ScopDomainSet.get_space();
  auto Ctx = ScopDomains.get_ctx();
  ScopSpace.dump();
  isl::space space(Ctx, 0, 3, 3);
  space = space.set_tuple_id(isl::dim::in, ScopDomainSet.get_tuple_id());
  space.dump();
  isl::map pwAffAsMap = isl::map::universe(space);
  pwAffAsMap = pwAffAsMap.project_out(isl::dim::out, 0, 1);
  //pwAffAsMap.dump();
  isl::local_space mapLocalSpace(pwAffAsMap.get_space());
  isl::constraint modulusConstraint = isl::constraint::alloc_equality(mapLocalSpace);
  modulusConstraint = modulusConstraint.set_coefficient_si(isl::dim::in, 1, 1);
  modulusConstraint = modulusConstraint.set_coefficient_si(isl::dim::out, 0, -2);
  modulusConstraint.dump();
  pwAffAsMap = pwAffAsMap.add_constraint(modulusConstraint);
  isl::constraint outputConstraint = isl::constraint::alloc_equality(mapLocalSpace);
  outputConstraint = outputConstraint.set_coefficient_si(isl::dim::in, 2, 1);
  outputConstraint = outputConstraint.set_coefficient_si(isl::dim::out, 1, -1);
  pwAffAsMap = pwAffAsMap.add_constraint(outputConstraint);
  pwAffAsMap = pwAffAsMap.project_out(isl::dim::out, 0, 1);
  pwAffAsMap.dump();
  isl::union_map pwAffAsUnionMap(pwAffAsMap);
  pwAffAsUnionMap.dump();
  auto aff = isl::multi_union_pw_aff::from_union_map(pwAffAsUnionMap);
  aff.dump();

  isl::map second = isl::map::universe(space);
  second = second.project_out(isl::dim::out, 0, 1);
  isl::local_space secondLocalSpace(second.get_space());
  isl::constraint secondModulus = isl::constraint::alloc_equality(secondLocalSpace);
  secondModulus = secondModulus.set_coefficient_si(isl::dim::in, 1, 1);
  secondModulus = secondModulus.set_coefficient_si(isl::dim::out, 0, -2);
  secondModulus = secondModulus.set_constant_si(1);
  second = second.add_constraint(secondModulus);
  isl::constraint secondOut = isl::constraint::alloc_equality(secondLocalSpace);
  secondOut = secondOut.set_coefficient_si(isl::dim::in, 2, 1);
  secondOut = secondOut.set_coefficient_si(isl::dim::out, 1, 1);
  second = second.add_constraint(secondOut);
  second = second.project_out(isl::dim::out, 0, 1);
  isl::union_map secondUnionMap(second);
  pwAffAsUnionMap = pwAffAsUnionMap.unite(secondUnionMap);
  pwAffAsUnionMap.dump();
  auto aff2 = isl::multi_union_pw_aff::from_union_map(pwAffAsUnionMap);
  aff2.dump();



//  isl::local_space localSpace(ScopSpace);
//  localSpace = localSpace.add_dims(isl::dim::div, 1);
//  localSpace.dump();
//  isl::aff aff(localSpace);
//  aff = aff.set_coefficient_si(isl::dim::in, 0, 3);
//  isl::aff aff2(localSpace);
//  aff2 = aff.set_coefficient_si(isl::dim::in, 0, 2);
//  aff = aff.add(aff2);
//  //aff = aff.floor();
//  aff.dump();
  //aff.dump();
}

// Define required info
void HMScheduleModifier::getAnalysisUsage(AnalysisUsage &AU) const {
  ScopPass::getAnalysisUsage(AU);
  AU.addRequired<polly::DependenceInfo>();
  AU.addPreserved<polly::DependenceInfo>();
  // This analysis blows up the schedule modifications of the former pass.
  // Hence, leave it away.
  //AU.addRequired<llvm::DependenceAnalysisWrapperPass>();
  //AU.addPreserved<llvm::DependenceAnalysisWrapperPass>();
}


Pass *polly::createHMScheduleModifierPass() { return new HMScheduleModifier(); }

INITIALIZE_PASS_BEGIN(HMScheduleModifier, "polly-hm-schedule-modifier",
                      "Polly - Custom pass by Hauke Mewes", false, false)
  INITIALIZE_PASS_DEPENDENCY(DependenceInfo)
  INITIALIZE_PASS_DEPENDENCY(ScopInfoRegionPass)
INITIALIZE_PASS_END(HMScheduleModifier, "polly-hm-schedule-modifier", "Polly - Custom pass by Hauke Mewes",
                    false, false)
#undef DEBUG_TYPE

