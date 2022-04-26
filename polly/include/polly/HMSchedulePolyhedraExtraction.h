//===------ HMSchedulePolyhedraExtraction.h ---------------------*- C++ -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// TODO!
//
//===----------------------------------------------------------------------===//
#ifndef POLLY_HM_SCHEDULE_POLYHEDRA_EXTRACTION_H
#define POLLY_HM_SCHEDULE_POLYHEDRA_EXTRACTION_H

#include "llvm/ADT/StringMap.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/Debug.h"

#include "isl/isl-noexceptions.h"

#include <random>

#define DEBUG_TYPE "hm-polyhedra-extraction"

/**
 * To understand the following, it is good to have read and somewhat understood
 * sections 3 and 4 of "Iterative Optimization in the Polyhedral Model:
 * Part I, One-Dimensional Time" by Pouchet et al.
 * I will base the following documentation on the code sample
 * of the matrix-vector multiplication that they provide, that is:
   for (i = 0; i <= n; i++) {
     S: s[i] = 0;
     for (j = 0; j <= n; j++) {
       R: s[i] = s[i] + a[i][j] * x[j];
     }
   }
 *
 * Furthermore, it is important to note that this code is adopted from
 * Polyite. https://github.com/stganser/polyite
 */


namespace polly {
namespace hmrtm {
/**
 * This class describes the characteristics of a single statement in a scop.
 * nrIt describes the number of domain iteration variables affecting this statement.
 * In the example case, there would be two StmtCoeffInfo objects, one for S and
 * one for R. As the statements are ordered in alphabetical order, R would appear
 * first in the list.
 * nrIt describes the number of domain variables that affect this statement,
 * itStart describes the offset of the first coefficient in the schedule row
 * for the domain variable coefficients. parStart does the same for the parameter
 * coefficients, and cstIdx for the constand coefficient. Lastly, the isl identifier
 * is stored.
 * In our example, this might look like this:
 * StmtCoeffInfo(identifier = R, itStart = 0, nrIt = 2, parStart = 3, cstIdx = 5)
 * StmtCoeffInfo(identifier = S, itStart = 2, nrIt = 1, parStart = 4, cstIdx = 6)
 */
class StmtCoeffInfo {
public:
  const int itStart;
  const int nrIt;
  const int parStart;
  const int cstIdx;
  const isl::id identifier;

  StmtCoeffInfo(const int itStart, const int nrIt, const int parStart,
                const int cstIdx, const isl::id &identifier);

};

/**
 * This class contains all the information about the schedule matrix rows.
 * It stores several different things. First of all, the
 * nrIt describes how many iteration multipliers exist in a program. For the
 * sample provided above, it would be three: 2 for the statement R, one for the
 * statement S. Then, nrParPS contains the number of external parameters of the
 * scop. For the example, this is one parameter, which is n.
 * The stmtInfo is explained above, and collects all the statements from this
 * scop. Dim describes the length of a schedule row, which is:
 * nrIt + nrParPS * nrStmts + 1 * nrStmts, resulting in 3 + 1 * 2 + 1 * 2 = 7.
 * The universe set is the base set for each schedule row with dim coefficients,
 * i.e. it is the set {[i0, i1, i2, i3, i4, i5, i6] : }
 */
class DomainCoeffInfo {
public:
  const int nrIt;
  const int nrParPS;
  const llvm::StringMap<StmtCoeffInfo> stmtInfo;
  const isl::set universe;
  const isl::union_set domain;
  const int nrStmts;
  const int dim;

  explicit DomainCoeffInfo(const isl::union_set &domain_);

private:

  int compNrIt(const isl::union_set &domain);

  int compNrParPS(const isl::union_set &domain);
  int compNrStmts(const isl::union_set &domain);

  llvm::StringMap<StmtCoeffInfo> compStmtInfo(const isl::union_set &domain);

  isl::set compUniverse(const isl::union_set &domain);

};
/**
 * This describes a dependence. It is stored in three different representations.
 * baseMap is the classical dependence relation of srcStmt -> targetStmt
 * weakConstr describes the constraint on the schedule coefficients that this
 * statement induces when it should be weakly satisfied in a schedule row.
 * strongConstr describes the constraint on the schedule coefficients that this
 * statement induces when it should be strongly satisfied in a schedule row.
 * Both constraints are represented as sets which contain all the values
 * that are still possible for the coefficients. This allows to combine
 * constraints by simply intersecting the sets. For example, to obtain the
 * one-dimensional schedule space that satisfies all dependencies strongly, one
 * only needs to intersect the univers set from DomainCoeffInfo with all
 * strongConstr of all dependencies.
 */
class Dependence {
public:
  // find a way to make these const.
  isl::basic_map baseMap;
  isl::basic_set weakConstr;
  isl::basic_set strongConstr;

  Dependence(const isl::basic_map &baseMap_, const isl::basic_set &weakConstr_,
             const isl::basic_set &strongConstr_);
  bool operator<(const Dependence& other) const;
};

/**
 * Calculate the dependencies and domain info from a scop, represented by the
 * following parameters.
 * @param params The params set. Can be obtained by domain.params() etc.
 * @param domain The domain of the scop.
 * @param sched The original schedule of the scop.
 * @param reads The reads of all statements.
 * @param writes The writes of all statements.
 * @return
 */
std::pair<DomainCoeffInfo,
                 llvm::SmallVector<Dependence,
                                   0>> calcDepsAndDomInfo(const isl::set &params,
                                                          isl::union_set domain,
                                                          isl::union_map sched,
                                                          const isl::union_map &reads,
                                                          const isl::union_map &writes);

llvm::SmallVector<isl::set, 0> constructMultiDimensionalSchedulePolyhedrons(const isl::set &params,
                                                                                   isl::union_set domain,
                                                                                   isl::union_map sched,
                                                                                   const isl::union_map &reads,
                                                                                   const isl::union_map &writes,
                                                                                   double prob2Carry,
                                                                                   std::pair<DomainCoeffInfo,
                                                                                             llvm::SmallVector<Dependence,
                                                                                                               0>>& depsAndDomInfo
);

/**
 * Construct a multi dimensional schedule space following the schedule space
 * creation described in: Iterative Schedule Optimization for Parallelization in
 * the Polyhedron Model by Ganser et al.
 * For the first parameters, see above.
 * @param params
 * @param domain
 * @param sched
 * @param reads
 * @param writes
 * @param prob2Carry Probability to carry at least one dependency in a schedule
 * dimension
 * @return a list containing sets describing the schedule dimensions.
 */
llvm::SmallVector<isl::set, 0> constructMultiDimensionalSchedulePolyhedrons(const isl::set &params,
                                                         isl::union_set domain,
                                                         isl::union_map sched,
                                                         const isl::union_map &reads,
                                                         const isl::union_map &writes,
                                                         double prob2Carry
    );
}
}
namespace llvm {
// move to cpp file, right now the linker complains that multiple
// definitions exist if this is not marked as static
llvm::raw_ostream& operator<<(llvm::raw_ostream &out,
                              const polly::hmrtm::StmtCoeffInfo &info);

llvm::raw_ostream& operator<<(llvm::raw_ostream &out,
                              const polly::hmrtm::DomainCoeffInfo &info);
}



namespace polly {
namespace hmrtm {
void testCalcDepsAndDomInfo(isl::ctx ctx);

}
}


// printing of dependence
#undef DEBUG_TYPE
#endif