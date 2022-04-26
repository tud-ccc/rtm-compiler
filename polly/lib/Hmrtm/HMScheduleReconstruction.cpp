//
// Created by hauke on 16.10.19.
//

#include "polly/HMScheduleReconstruction.h"

using namespace polly::hmrtm;
isl::constraint polly::hmrtm::createScheduleConstraint(const std::vector<int> &ScheduleVec,
                                                       const polly::hmrtm::StmtCoeffInfo &Stmt,
                                                       isl::local_space &LocalSpace,
                                                       int Dimension) {
  auto Ctx = LocalSpace.get_ctx();
  auto Constraint = isl::constraint::alloc_equality(LocalSpace);
  Constraint = Constraint.set_coefficient_val(
      isl::dim::out, Dimension, isl::val::negone(Ctx));
  // ok, do that the next time. This is too hard for my brain right now.
  for(int i = 0; i < Stmt.nrIt; ++i) {
    Constraint = Constraint.set_coefficient_si(isl::dim::in, i, ScheduleVec[i + Stmt.itStart]);
  }
  int nrParams = LocalSpace.dim(isl::dim::param);
  for(int i = 0; i < nrParams; ++i) {
    Constraint = Constraint.set_coefficient_si(isl::dim::param, i, ScheduleVec[i + nrParams]);
  }
  Constraint = Constraint.set_constant_si(ScheduleVec[Stmt.cstIdx]);
  return Constraint;
}
isl::local_space polly::hmrtm::createLocalSpace(isl::ctx Ctx,
                                                const polly::hmrtm::DomainCoeffInfo &DomInfo,
                                                const polly::hmrtm::StmtCoeffInfo &stmt,
                                                const unsigned scheduleDimensions) {
  auto lSp = isl::local_space(isl::space(
      Ctx,
      static_cast<unsigned>(DomInfo.nrParPS),
      static_cast<unsigned>(stmt.nrIt),
      scheduleDimensions));
  lSp = lSp.set_tuple_id(isl::dim::in, stmt.identifier);
  auto params = DomInfo.domain.params();
  for(int i = 0; i < DomInfo.nrParPS; ++i) {
    lSp = lSp.set_dim_id(isl::dim::param, i, params.get_dim_id(isl::dim::param, i));
  }
  return lSp;
}
isl::union_map polly::hmrtm::createMultiDimensionalSchedule(isl::ctx Ctx,
                                                            const polly::hmrtm::DomainCoeffInfo &DomInfo,
                                                            const std::vector<
                                                                std::vector<int>> &ScheduleVecs) {
  auto Dims = ScheduleVecs.size();
  std::vector<isl::map> mapList;
  for(auto &entry: DomInfo.stmtInfo) {
    auto &stmt = entry.second;
    auto localSpace = createLocalSpace(Ctx, DomInfo, stmt, Dims);
    auto ScheduleMap = isl::map::universe(localSpace.get_space());
    for(decltype(Dims) i = 0; i < Dims; ++i) {
      // TODO: uncomment
      auto constraint = createScheduleConstraint(ScheduleVecs[i], stmt, localSpace,
                                                 static_cast<int>(i));
      ScheduleMap = ScheduleMap.add_constraint(constraint);
    }
    mapList.push_back(ScheduleMap);
  }
  auto ScheduleMap = isl::union_map(mapList[mapList.size() - 1]);
  mapList.pop_back();
  for(auto &nextStmtMap: mapList) {
    ScheduleMap = ScheduleMap.unite(nextStmtMap);
  }
  return ScheduleMap;
}
