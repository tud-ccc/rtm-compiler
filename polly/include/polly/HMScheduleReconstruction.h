//===------ HMScheduleReconstruction.h --------------------------*- C++ -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This allows to convert schedule matrices to isl::union_maps which can be
// used as schedules in polly.
//
//===----------------------------------------------------------------------===//

#ifndef POLLY_HM_SCHEDULE_RECONSTRUCTION_H
#define POLLY_HM_SCHEDULE_RECONSTRUCTION_H

#include "polly/HMSchedulePolyhedraExtraction.h"

namespace polly {
namespace hmrtm {

isl::constraint createScheduleConstraint(
    const std::vector<int> &ScheduleVec,
    const StmtCoeffInfo &Stmt,
    isl::local_space &LocalSpace,
    int Dimension);

isl::local_space createLocalSpace(isl::ctx Ctx,
    const DomainCoeffInfo &DomInfo,
    const StmtCoeffInfo &stmt,
    const unsigned scheduleDimensions);

/**
 * This converts a schedule matrix and the corresponding domain coefficient info
 * back to an isl::union_map
 * @param Ctx isl context.
 * @param DomInfo
 * @param ScheduleVecs
 * @return schedule as union_map, to be used in polly.
 */
isl::union_map createMultiDimensionalSchedule(isl::ctx Ctx,
    const DomainCoeffInfo &DomInfo,
    const std::vector<std::vector<int>> &ScheduleVecs);

}
}

#endif //POLLY_HM_SCHEDULE_RECONSTRUCTION_H
