#include "polly/HMSchedulePolyhedraExtraction.h"
#include "polly/HMScheduleShiftCalculator.h"
#include "polly/HMScheduleReconstruction.h"
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

#define DEBUG_TYPE "polly-opt-hm-schedule-search"


using namespace llvm;
using namespace polly;

static cl::opt<bool>
    HMPollyEnableFullSearch("polly-hm-enable-full-search", cl::desc("Enable the full shift search."),
                 cl::init(false), cl::ZeroOrMore);

template<typename vec>
using VecType = std::vector<vec>;

namespace {
class HMRacetrackScheduleSpaceSearcher : public ScopPass {
public:
  static char ID;

  explicit HMRacetrackScheduleSpaceSearcher() : ScopPass(ID) {}

  /// Optimize the schedule of the SCoP @p S for SPM utilizing racetrack memory.
  bool runOnScop(Scop &S) override;

  /// Register all analyses and transformations required
  void getAnalysisUsage(AnalysisUsage &AU) const override;
};
} // namespace

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

char HMRacetrackScheduleSpaceSearcher::ID = 0;

static isl::set assignValueToSet(isl::set set, isl::dim type, int pos, int value) {
  isl::local_space ConstraintSpace(set.get_space());
  auto Constraint = isl::constraint::alloc_equality(ConstraintSpace);
  Constraint = Constraint.set_constant_si(value);
  Constraint = Constraint.set_coefficient_si(type, pos, -1);
  return set.add_constraint(Constraint);
}

template<typename T>
static void computeScheduleCombinations(polly::hmrtm::DomainCoeffInfo &DomInfo,
    const VecType<isl::set> &FourierMotzkinSets,
    const T& ValuesToTry,
    VecType<VecType<int>> &Results,
    VecType<int> CurrentSequence, int newValue) {
  // assign newValue as first parameter of fourier-motzkin set
  isl::set FMTestSet = FourierMotzkinSets[CurrentSequence.size()];
  FMTestSet = assignValueToSet(FMTestSet, isl::dim::set, 0, newValue);
  if(!FMTestSet.is_empty()) {
    CurrentSequence.push_back(newValue);
    if(static_cast<int>(CurrentSequence.size()) == DomInfo.dim) {
      Results.emplace_back(std::move(CurrentSequence));
    } else {
      for(auto V: ValuesToTry) {
        computeScheduleCombinations(DomInfo, FourierMotzkinSets, ValuesToTry, Results, CurrentSequence, V);
      }
    }
  }
}

static void calculateShiftsForAllPermutations(polly::hmrtm::DomainCoeffInfo &DomInfo,
    Scop &S, VecType<VecType<VecType<int>>> &PossibleCoefficients,
   const VecType<VecType<int>> &ScheduleVecs) {
  if(ScheduleVecs.size() == PossibleCoefficients.size()) {
    auto Schedule = polly::hmrtm::createMultiDimensionalSchedule(S.getIslCtx(), DomInfo, ScheduleVecs);
    auto shifts = polly::hmrtm::calculateNumberOfShifts(S, Schedule);
    dbgs() << "Shifts: " << shifts << "\n";
  } else {
    auto size = ScheduleVecs.size();
    for(auto &vec: PossibleCoefficients[size]) {
      VecType<VecType<int>> ScheduleVecsCopy(ScheduleVecs);
      ScheduleVecsCopy.push_back(vec);
      calculateShiftsForAllPermutations(DomInfo, S, PossibleCoefficients, ScheduleVecsCopy);
    }
  }
}

template<typename T>
VecType<VecType<isl::set>> scheduleSpaceFourierMotzkinProjection(const T& ScheduleSpace) {
  auto Size = ScheduleSpace.size();
  VecType<VecType<isl::set>> FMProjection(ScheduleSpace.size());
  for(decltype(Size) i = 0; i < Size; ++i) {
    isl::set LastSet = ScheduleSpace[i];
    while(LastSet.dim(isl::dim::set) > 0) {
      FMProjection[i].push_back(LastSet);
      LastSet = LastSet.project_out(isl::dim::set, 0, 1);
    }
  }
  return FMProjection;
}

template<typename P>
void createRandomScheduleRowHelper(const P &PossibleValueGenerator, VecType<isl::set> &FMRow,
    VecType<int> &CurrentSchedule) {
  auto Pos = CurrentSchedule.size();
  // Schedule Row is complete
  if(Pos == FMRow.size()) {
    return;
  }
  int nextInt = PossibleValueGenerator();
  auto fmAssigned = assignValueToSet(FMRow[Pos], isl::dim::set, 0, nextInt);
  if(!fmAssigned.is_empty()) {
    CurrentSchedule.push_back(nextInt);
    createRandomScheduleRowHelper(PossibleValueGenerator, FMRow, CurrentSchedule);
  }
}

template<typename P>
VecType<int> createRandomScheduleRow(const P &PossibleValueGenerator, VecType<isl::set> &FMRow) {
  VecType<int> Result;
  while(Result.size() < FMRow.size()) {
    Result.clear();
    createRandomScheduleRowHelper(PossibleValueGenerator, FMRow, Result);
  }
  return Result;
}

template<typename P>
VecType<VecType<int>> createRandomSchedule(const P& PossibleValues, VecType<VecType<isl::set>>& FMProjection) {

  //initialize rng
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> SizeSelect(0, PossibleValues.size() - 1);
  auto PossibleValueGenerator = [&SizeSelect, &PossibleValues, &rng]() {
    return PossibleValues[SizeSelect(rng)];
  };

  VecType<VecType<int>> RandomSchedule(FMProjection.size());
  for(auto &FMRow: FMProjection) {
    RandomSchedule.emplace_back(createRandomScheduleRow(PossibleValueGenerator, FMRow));
  }
  return RandomSchedule;
}

bool HMRacetrackScheduleSpaceSearcher::runOnScop(Scop &S) {

  LLVM_DEBUG(dbgs() << "In Racetrack Schedule Optimizer!\n");
  auto shifts = polly::hmrtm::calculateNumberOfShifts(S, S.getSchedule());
  LLVM_DEBUG(dbgs() << "Shifts: " << shifts << "\n");
//  LLVM_DEBUG({
//    dbgs() << "Program calculated for shifts";
//    auto P = polly::hmrtm::createShiftProgram(S);
//    for(auto &S: P.Statements) {
//      S.Domain.dump();
//      for(auto &A: S.MemoryAccesses) {
//        A.dump();
//      }
//    }
//  });

  auto Domains = S.getDomains();
  auto Params = Domains.params();
  auto Schedule = S.getSchedule();
  auto Reads = S.getReads().gist_domain(Domains);
  auto Writes = S.getWrites().gist_domain(Domains);
//  LLVM_DEBUG({
//    dbgs() << "Debuggin params to schedule polyhedron construction\n";
//    Domains.dump();
//    Params.dump();
//    Schedule.dump();
//    S.getReads().dump();
//    S.getWrites().dump();
//    auto DepAndDomInfo = polly::hmrtm::calcDepsAndDomInfo(
//        Params, Domains, Schedule, Reads, Writes);
//    auto schedulePoly = DepAndDomInfo.first.universe;
//    dbgs() << "Deps: " << DepAndDomInfo.second.size() << "\n";
//    dbgs() << "Dumping deps...\n";
//    for(auto &dep: DepAndDomInfo.second) {
//      //dep.strongConstr.dump();
//      dbgs() << dep.strongConstr.to_str() << ",,,,";
//      schedulePoly = schedulePoly.intersect(dep.strongConstr);
//    }
//    dbgs() << "\nDumping schedule poly...\n";
//    schedulePoly.dump();
//  });
  auto DepAndDomInfo = polly::hmrtm::calcDepsAndDomInfo(
      Params, Domains, Schedule, Reads, Writes);
  auto ScheduleSpace = polly::hmrtm::constructMultiDimensionalSchedulePolyhedrons(
      Params, Domains, Schedule, Reads, Writes, 1.0, DepAndDomInfo);


//  {
//    auto FMProjection = scheduleSpaceFourierMotzkinProjection(ScheduleSpace);
//    for(auto &V: FMProjection) {
//      for(auto &Proj: V) {
//        Proj.dump();
//      }
//    }
//    llvm::dbgs() << "Starting random schedule generation...";
//    VecType<int> PossibleValues{-1, 0, 1};
//    auto Sched = createRandomSchedule(PossibleValues, FMProjection);
//    auto ScheduleMap = polly::hmrtm::createMultiDimensionalSchedule(S.getIslCtx(),
//        DepAndDomInfo.first, Sched);
//    shifts = polly::hmrtm::calculateNumberOfShifts(S, ScheduleMap);
//    llvm::dbgs() << "New shifts: " << shifts << "\n";
//  }
  auto newSchedule = isl::union_map(S.getIslCtx(), "{ Stmt2[i0, i1, i2] -> [82 + 2i0, -68 + 5i0 + i1 + i2, 97 + 2i0 + 2i1 + 2i2, -86 + 4i0 - 4i1]; Stmt5[i0, i1, i2] -> [83 + 2i0, -65 + 5i0 + i1 + i2, 99 + 2i0 - i1 + 2i2, -84 + 4i0 + 5i1 + 3i2] }");
  std::vector<isl::id> ids;
  auto MapList = Schedule.get_map_list();
  for(int i = 0; i < MapList.size(); ++i) {
    auto Map = MapList.get_at(i);
    auto id = Map.get_tuple_id(isl::dim::in);
    ids.push_back(id);
  }

  MapList = newSchedule.get_map_list();
  newSchedule = isl::union_map::empty(Schedule.get_space());
  for(int i = 0; i < MapList.size(); ++i) {
    auto Map = MapList.get_at(i);
    auto id = Map.get_tuple_name(isl::dim::in);
    for(auto &NewId: ids) {
      if(NewId.get_name() == id) {
        LLVM_DEBUG(llvm::dbgs() << "Exchanged id! " << NewId.to_str() << "\n");
        Map = Map.set_tuple_id(isl::dim::in, NewId);
      }
    }
    newSchedule = newSchedule.add_map(Map);
  }
  llvm::dbgs() << "Wow";
  newSchedule.intersect_domain(Domains).dump();
  S.setSchedule(newSchedule);
  printAstForScop(S);


  if(HMPollyEnableFullSearch) {

    for(auto &d: ScheduleSpace) {
      d.dump();
    }
    std::array<int, 3> PossibleIndices = {-1, 0, 1};
    VecType<VecType<isl::set>> FourierMotzkinSets(ScheduleSpace.size());
    for(decltype(ScheduleSpace.size()) i = 0, size = ScheduleSpace.size(); i < size; ++i) {
      isl::set LastSet = ScheduleSpace[i];
      while(LastSet.dim(isl::dim::set) > 0) {
        FourierMotzkinSets[i].push_back(LastSet);
        LastSet = LastSet.project_out(isl::dim::set, 0, 1);
      }
    }
//    for(auto &Projected: FourierMotzkinSets[0]) {
//      Projected.dump();
//    }
    // First interesting experiment:
    // Do some tests on all combinations:
    auto &DomInfo = DepAndDomInfo.first;
    VecType<VecType<VecType<int>>> PossibleCoefficientsPerDim(ScheduleSpace.size());
    for(decltype(ScheduleSpace.size()) ScheduleDimension = 0;
    ScheduleDimension < ScheduleSpace.size();
    ++ScheduleDimension) {
      for(auto V: PossibleIndices) {
        computeScheduleCombinations(DomInfo, FourierMotzkinSets[ScheduleDimension], PossibleIndices,
            PossibleCoefficientsPerDim[ScheduleDimension], VecType<int>(), V);
      }
    }
    LLVM_DEBUG({
      dbgs() << "Results for first dimension level.\n";
      calculateShiftsForAllPermutations(DomInfo,S, PossibleCoefficientsPerDim, VecType<VecType<int>>{});
      dbgs() << "Results end.\n";
    });


  }
  
  return false;
}

void HMRacetrackScheduleSpaceSearcher::getAnalysisUsage(AnalysisUsage &AU) const {
  ScopPass::getAnalysisUsage(AU);
  AU.addRequired<DependenceInfo>();
  AU.addPreserved<DependenceInfo>();
}

Pass *polly::createHMRacetrackScheduleSpaceSearcherPass() {
  return new HMRacetrackScheduleSpaceSearcher();
}

INITIALIZE_PASS_BEGIN(HMRacetrackScheduleSpaceSearcher,
    "polly-opt-hm-schedule-search",
    "Polly Racetrack Schedule Space Searcher - Custom pass by Hauke Mewes", false, false)
INITIALIZE_PASS_DEPENDENCY(DependenceInfo)
INITIALIZE_PASS_DEPENDENCY(ScopInfoRegionPass)
INITIALIZE_PASS_END(HMRacetrackScheduleSpaceSearcher,
    "polly-opt-hm-schedule-search",
    "Polly Racetrack Schedule Space Searcher - Custom pass by Hauke Mewes", false, false)

#undef DEBUG_TYPE