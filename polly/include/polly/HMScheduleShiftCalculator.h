//===------ HMScheduleShiftCalculator.h -------------------------*- C++ -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file provides a function to calculate the shifts on a racetrack memory
// model with 64 entries in each bench, assuming that the entire memory fits
// into the rtm.
// TODO: make the model more adaptive.
//
//===----------------------------------------------------------------------===//
#ifndef POLLY_HMSCHEDULESHIFTCALCULATOR_H
#define POLLY_HMSCHEDULESHIFTCALCULATOR_H

#include "polly/ScopInfo.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/SmallSet.h"
#include "llvm/ADT/StringSet.h"
#include "llvm/ADT/StringMap.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/raw_ostream.h"

#include <vector>
#include <regex>
#include <tuple>

#define DEBUG_TYPE "polly-hm-trace-generation"

template<typename vec>
using VecType = std::vector<vec>;

namespace polly {
namespace hmrtm {

struct ShiftConfig {
  std::size_t RTMSize;
  bool TraceFile;
  std::string KernelName;
  std::string OutputDirectory;
  std::string Position;

  ShiftConfig(std::size_t RTMSize_, bool TraceFile_) :
  RTMSize(RTMSize_), TraceFile(TraceFile_), KernelName(""), OutputDirectory(""),
  Position("") {}

  explicit ShiftConfig(std::size_t RTMSize_) : ShiftConfig(RTMSize_, false) {

  }
  ShiftConfig() : ShiftConfig(64, false) {

  }

  ShiftConfig& traceFile(bool TraceFile_) {
    TraceFile = TraceFile_;
    return *this;
  }

  ShiftConfig& rtmSize(std::size_t Size) {
    RTMSize = Size;
    return *this;
  }

  ShiftConfig& kernelName(const std::string& KernelName_) {
    KernelName = KernelName_;
    return *this;
  }
  ShiftConfig& outputDirectory(const std::string& OutputDirectory_) {
    OutputDirectory = OutputDirectory_;
    return *this;
  }
  ShiftConfig& position(const std::string &Position_) {
    Position = Position_;
    return *this;
  }

  template<typename Stream>
  void print(Stream &Out) const {
    Out << "RTMSize = " << RTMSize
    << "\nTraceFile = " << TraceFile
    << "\nKernelName = " << KernelName
    << "\nOutputDirectory = " << OutputDirectory
    << "\nPosition = " << Position << "\n";
  }
};

struct Statement {
  isl::set Domain;
  VecType<std::pair<isl::map, bool>> MemoryAccesses;

  explicit Statement(const isl::set &Domain);

  Statement(Statement &&s) = default;
  Statement(const Statement &s);

  void addAccess(const isl::map &Access, bool IsRead);
};

struct ShiftCalculationProgram {
  isl::union_set Domains;
  isl::set Context;
  //isl::union_map Schedule;
  VecType<Statement> Statements;
  llvm::SmallSet<std::string, 5> Variables;

  explicit ShiftCalculationProgram(const isl::union_set &Domains, const isl::set &Context);

  void addStatement(Statement &&Stmt);
};

ShiftCalculationProgram createShiftProgram(polly::Scop &S);


std::string astExprToStr(const isl::ast_expr &ast);

std::string astNodeToStr(const isl::ast_node &ast);

/**
 * This outputs the preliminaries to the C++ file that is used to calculates the
 * shifts. This includes some elaborate C++ template magic to create a class
 * with an arbitrary number of method arguments.
 * @tparam Stream Any type that supports the << operator on char*.
 * @param Out An output stream.
 */
template<typename Stream>
void printPreliminariesToStream(Stream &Out) {
//  Out << "#include <utility>\n"
//         "#include <iostream>\n"
//         "#include <tuple>\n"
//         "#include <vector>\n"
//         "#include <algorithm>\n"
//         "using namespace std;\n"
//         "\n"
//         "\n"
//         "static int floord(int n, int d) {return ((n < 0) ? (n - d + 1) : n) / d;}\n"
//         "template<typename Tuple, std::size_t index>\n"
//         "struct MultiplyTuple {\n"
//         "  static std::size_t multiply(const Tuple &tuple) {\n"
//         "    return std::get<std::tuple_size<Tuple>::value  - index>(tuple) * MultiplyTuple<Tuple, index - 1>::multiply(tuple);\n"
//         "  }\n"
//         "};\n"
//         "\n"
//         "// base case has special value\n"
//         "template<typename Tuple>\n"
//         "struct MultiplyTuple<Tuple, 1> {\n"
//         "  static std::size_t multiply(const Tuple &tuple) {\n"
//         "    return 1;\n"
//         "  }\n"
//         "};\n"
//         "\n"
//         "template<typename Tuple>\n"
//         "static std::size_t multiplyAllButOne(const Tuple &tuple) {\n"
//         "  return MultiplyTuple<Tuple, std::tuple_size<Tuple>::value>::multiply(tuple);\n"
//         "}\n"
//         "\n"
//         "template<typename TRM, TRM rtmSize, typename... Types>\n"
//         "struct TraceKeeper {\n"
//         "  std::tuple<Types...> sizes;\n"
//         "  std::size_t innermostDimensionSize;\n"
//         "  std::vector<std::size_t> rtmHeads;\n"
//         "  std::size_t shifts;\n"
//         "\n"
//         "  //std::tuple<int, int> test;\n"
//         "\n"
//         "  explicit TraceKeeper(std::tuple<Types...> sizes_)\n"
//         "  : sizes(std::move(sizes_)),\n"
//         "  innermostDimensionSize(std::get<std::tuple_size<decltype(sizes)>::value - 1>(sizes) / rtmSize + 1),\n"
//         "  rtmHeads(multiplyAllButOne(sizes) * innermostDimensionSize),\n"
//         "  shifts(0) {\n"
//         "    //std::cout << \"Total size: \" << rtmHeads.size() << \"\\n\";\n"
//         "  }\n"
//         "\n"
//         "  void access(Types... values) {\n"
//         "    implAccess<0>(0, values...);\n"
//         "  }\n"
//         "\n"
//         "private:\n"
//         "  template<std::size_t depth, typename ...T>\n"
//         "  void implAccess(std::size_t currentIndex, int a, T... vals) {\n"
//         "    //std::cout << a << \" \";\n"
//         "    auto currentMult = std::get<depth>(sizes);\n"
//         "    implAccess<depth + 1>(a + currentIndex * currentMult, vals...);\n"
//         "  }\n"
//         "  template<std::size_t depth>\n"
//         "  void implAccess(std::size_t currentIndex, int a) {\n"
//         "    // depth should point to the last tuple element index\n"
//         "    static_assert(depth == std::tuple_size<decltype(sizes)>::value - 1, \"Recursion depth must equal to the number of tuple elements.\");\n"
//         "    // This is the base case, a is the outer most index now\n"
//         "    // so truncate by rtm memory size\n"
//         "    auto currentMult = innermostDimensionSize;\n"
//         "    auto finalIndex = static_cast<std::size_t>(a / rtmSize) + currentMult * currentIndex;\n"
//         "    auto oldHeadPosition = rtmHeads[finalIndex];\n"
//         "    auto newHeadPosition = static_cast<TRM>(a % rtmSize);\n"
//         "    auto addedShifts = std::abs(static_cast<int>(newHeadPosition - oldHeadPosition));\n"
//         "    shifts += addedShifts;\n"
//         "    rtmHeads[finalIndex] = newHeadPosition;\n"
//         "    //std::cout << a << \":\" << addedShifts << \";\" << finalIndex << \"\\n\";\n"
//         "  }\n"
//         "\n"
//         "};";
  Out << "#include <utility>\n"
         "#include <iostream>\n"
         "#include <fstream>\n"
         "#include <tuple>\n"
         "#include <vector>\n"
         "#include <algorithm>\n"
         "#include <iomanip>\n"
         "using namespace std;\n"
         "\n"
         "\n"
         "static int floord(int n, int d) {return ((n < 0) ? (n - d + 1) : n) / d;}\n"
         "template<typename Tuple, std::size_t index>\n"
         "struct MultiplyTuple {\n"
         "  static std::size_t multiply(const Tuple &tuple) {\n"
         "    return std::get<std::tuple_size<Tuple>::value  - index>(tuple) * MultiplyTuple<Tuple, index - 1>::multiply(tuple);\n"
         "  }\n"
         "};\n"
         "\n"
         "// base case has special value\n"
         "template<typename Tuple>\n"
         "struct MultiplyTuple<Tuple, 1> {\n"
         "  static std::size_t multiply(const Tuple &tuple) {\n"
         "    return 1;\n"
         "  }\n"
         "};\n"
         "\n"
         "template<typename Tuple>\n"
         "static std::size_t multiplyAllButOne(const Tuple &tuple) {\n"
         "  return MultiplyTuple<Tuple, std::tuple_size<Tuple>::value>::multiply(tuple);\n"
         "}\n"
         "\n"
         "using RowCol = std::pair<std::size_t, std::size_t>;\n"
         "\n"
         "template<typename TRM, TRM rtmSize, typename... Types>\n"
         "struct TraceKeeper {\n"
         "  std::tuple<Types...> sizes;\n"
         "  std::size_t innermostDimensionSize;\n"
         "  std::vector<std::size_t> rtmHeads;\n"
         "  std::size_t shifts;\n"
         "\n"
         "  //std::tuple<int, int> test;\n"
         "\n"
         "  explicit TraceKeeper(std::tuple<Types...> sizes_)\n"
         "      : sizes(std::move(sizes_)),\n"
         "        innermostDimensionSize((std::get<std::tuple_size<decltype(sizes)>::value - 1>(sizes) + rtmSize - 1) / rtmSize),\n"
         "        rtmHeads(multiplyAllButOne(sizes) * innermostDimensionSize),\n"
         "        shifts(0) {\n"
         "    //std::cout << \"Total size: \" << rtmHeads.size() << \"\\n\";\n"
         "  }\n"
         "\n"
         "  void access(Types... values) {\n"
         "    //implAccess<0>(0, values...);\n"
         "    access(computeRowCol<0>(0, values...));\n"
         "\n"
         "  }\n"
         "\n"
         "  template<typename Stream>\n"
         "  void print(Stream &out) {\n"
         "    out << \"  innermostDimensionSize: \" << innermostDimensionSize\n"
         "    << \"\\n  rtmHeads.size(): \" << rtmHeads.size()\n"
         "    << \"\\n  shifts: \" << shifts << \"\\n\";\n"
         "  }\n"
         "\n"
         "protected:\n"
         "\n"
         "  void access(const RowCol &rowCol) {\n"
         "    std::size_t row;\n"
         "    std::size_t col;\n"
         "    std::tie(row, col) = rowCol;\n"
         "    auto oldHeadPosition = rtmHeads[row];\n"
         "    auto newHeadPosition = static_cast<TRM>(col);\n"
         "    auto addedShifts = std::abs(static_cast<int>(newHeadPosition - oldHeadPosition));\n"
         "    shifts += addedShifts;\n"
         "    rtmHeads[row] = newHeadPosition;\n"
         "  }\n"
         "\n"
         "  template<std::size_t depth, typename ...T>\n"
         "  RowCol\n"
         "  computeRowCol(std::size_t currentIndex, int a, T... vals) {\n"
         "    auto currentMult = std::get<depth>(sizes);\n"
         "    return computeRowCol<depth + 1>(a + currentIndex * currentMult, vals...);\n"
         "  }\n"
         "\n"
         "  template<std::size_t depth>\n"
         "  RowCol\n"
         "  computeRowCol(std::size_t currentIndex, int a) {\n"
         "    static_assert(depth == std::tuple_size<decltype(sizes)>::value - 1, \"Recursion depth must equal to the number of tuple elements.\");\n"
         "    auto currentMult = innermostDimensionSize;\n"
         "    auto row = static_cast<std::size_t>(a / rtmSize) + currentMult * currentIndex;\n"
         "    auto col = static_cast<std::size_t>(a % rtmSize);\n"
         "    return std::make_pair(row, col);\n"
         "  }\n"
         "\n"
         "//  template<std::size_t depth, typename ...T>\n"
         "//  void implAccess(std::size_t currentIndex, int a, T... vals) {\n"
         "//    //std::cout << a << \" \";\n"
         "//    auto currentMult = std::get<depth>(sizes);\n"
         "//    implAccess<depth + 1>(a + currentIndex * currentMult, vals...);\n"
         "//  }\n"
         "//  template<std::size_t depth>\n"
         "//  void implAccess(std::size_t currentIndex, int a) {\n"
         "//    // depth should point to the last tuple element index\n"
         "//    static_assert(depth == std::tuple_size<decltype(sizes)>::value - 1, \"Recursion depth must equal to the number of tuple elements.\");\n"
         "//    // This is the base case, a is the outer most index now\n"
         "//    // so truncate by rtm memory size\n"
         "//    auto currentMult = innermostDimensionSize;\n"
         "//    auto finalIndex = static_cast<std::size_t>(a / rtmSize) + currentMult * currentIndex;\n"
         "//    auto oldHeadPosition = rtmHeads[finalIndex];\n"
         "//    auto newHeadPosition = static_cast<TRM>(a % rtmSize);\n"
         "//    auto addedShifts = std::abs(static_cast<int>(newHeadPosition - oldHeadPosition));\n"
         "//    shifts += addedShifts;\n"
         "//    rtmHeads[finalIndex] = newHeadPosition;\n"
         "//    //std::cout << a << \":\" << addedShifts << \";\" << finalIndex << \"\\n\";\n"
         "//  }\n"
         "};\n"
         "\n"
         "\n"
         "enum AccessType {\n"
         "  READ, WRITE\n"
         "};\n"
         "\n"
         "template<std::size_t rtmSize, typename TracePrinterObj, typename... Types>\n"
         "struct TracePrinter : TraceKeeper<std::size_t, rtmSize, Types...> {\n"
         "private:\n"
         "  std::ofstream traceFile;\n"
         "  TracePrinterObj* printer;\n"
         "public:\n"
         "\n"
         "  std::size_t memOffset;\n"
         "  using Base = TraceKeeper<std::size_t, rtmSize, Types...>;\n"
         "  TracePrinter(std::tuple<Types...> sizes_,\n"
         "      std::size_t memOffset_,\n"
         "      TracePrinterObj* printer_) : Base(sizes_), traceFile(), printer(printer_), memOffset(memOffset_) {\n"
         "  }\n"
         "\n"
         "  std::size_t memRange() const {\n"
         "    return Base::rtmHeads.size() * rtmSize;\n"
         "  }\n"
         "\n"
         "  std::size_t addressAfter() const {\n"
         "    return memRange() + memOffset;\n"
         "  }\n"
         "\n"
         "  void accessAndPrint(AccessType access, Types... values) {\n"
         "    auto rowCol = Base::template computeRowCol<0>(0, values...);\n"
         "    Base::access(rowCol);\n"
         "    printer->printToStream(rowCol, access, memOffset);\n"
         "    //TracePrinterObj::printToStream(std::cout, rowCol, memOffset);\n"
         "  }\n"
         "\n"
         "  template<typename Stream>\n"
         "  void print(Stream &out) {\n"
         "    out << \"  innermostDimensionSize: \" << Base::innermostDimensionSize\n"
         "        << \"\\n  rtmHeads.size(): \" << Base::rtmHeads.size()\n"
         "        << \"\\n  shifts: \" << Base::shifts\n"
         "        << \"\\n  memOffset: \" << memOffset\n"
         "        << \"\\n  memRange: \" << memRange() << \"\\n\";\n"
         "  }\n"
         "\n"
         "};\n"
         "\n"
         "//template<typename StreamType>\n"
         "struct TracePrinterRTSim {\n"
         "  const std::size_t rtmSize;\n"
         "  std::size_t counter;\n"
         "  std::ofstream traceFile;\n"
         "\n"
         "  explicit TracePrinterRTSim(std::size_t rtmSize_, const std::string &traceFileName) :rtmSize(rtmSize_), counter(0), traceFile() {\n"
         "    traceFile.open(traceFileName);\n"
         "    traceFile << \"\\n\";\n"
         "  }\n"
         "\n"
         "  void printToStream(const RowCol &rowCol, const AccessType access, const std::size_t memOffset) {\n"
         "    std::size_t row;\n"
         "    std::size_t col;\n"
         "    std::tie(row, col) = rowCol;\n"
         "    auto address = row * rtmSize + col + memOffset;\n"
         "    address = 2 * address * 32;\n"
         "    ++counter;\n"
         "    auto accessStr = access == AccessType::READ ? \"R\" : \"W\";\n"
         "    traceFile << std::dec << counter << \" \" << accessStr << \" 0x\" << std::hex << address\n"
         "    << \" 00e8da000000e8d5eb05004883c408c331ed4989d15e4889e24883e4f0505449c7c0a00a400048c7c1e00a400048c7c7a4024000e8a7010000f490904883ec08 0\\n\";\n"
         "  }\n"
         "\n"
         "  ~TracePrinterRTSim() {\n"
         "    traceFile.close();\n"
         "  }\n"
         "};";
}

/**
 * This completes the shift calculation code.
 * @tparam Stream
 * @param Out
 * @param Program
 * @param Schedule
 */
template<typename Stream>
static void printProgramToStream(Stream &Out, ShiftCalculationProgram &Program, const isl::union_map &Schedule,
    const ShiftConfig &Config, const std::string &TraceFile) {
  bool CreateTraceFile = Config.TraceFile;
  printPreliminariesToStream(Out);
  Out << "struct Runner {\n";
  if(CreateTraceFile) {
    Out << "  TracePrinterRTSim printer_;\n";
  }
  llvm::StringSet<> VarNames;
  for(auto &Stmt: Program.Statements) {
    for(auto &Access: Stmt.MemoryAccesses) {
      int dims = Access.first.dim(isl::dim::out);
      // filter out phi nodes
      if(dims == 0) {
        continue;
      }
      std::string VarName = Access.first.get_tuple_name(isl::dim::out);
      // This checks for non-existence of VarName
      if(VarNames.find(VarName) == VarNames.end()) {
        if(!CreateTraceFile) {
          Out << "  TraceKeeper<int, " << Config.RTMSize <<", ";
        } else {
          Out << "  TracePrinter<" << Config.RTMSize << ", TracePrinterRTSim, ";
        }
        for(int i = 0; i < dims; ++i) {
          Out << "int";
          if(i < dims - 1) {
            Out << ", ";
          }
        }
        Out << "> "
          << VarName << ";\n";
        VarNames.insert(VarName);
      }
    }
  }

  // next, find max of each var name
  //Out << "Before max finder!\n";
  VecType<std::string> MemRefNames(VarNames.size());
  VecType<VecType<long>> MemRefValues(VarNames.size());
  {
    decltype(MemRefNames.size()) i = 0;
    for(auto &MemRefName: VarNames) {
      MemRefNames[i] = MemRefName.first().str();
      isl::set a;
      bool found = false;
      for(auto &Stmt: Program.Statements) {
        for(auto &Access: Stmt.MemoryAccesses) {
          auto MemAccessRefName = Access.first.get_tuple_name(isl::dim::out);
          if(MemRefName.first() == MemAccessRefName) {
            if(!found) {
              found = true;
              a = Stmt.Domain.apply(Access.first);
            } else {
              a = a.unite(Stmt.Domain.apply(Access.first));
            }
          }
        }
      }
      auto dimensions = a.dim(isl::dim::set);
      //a.dump();
      for(decltype(dimensions) j = 0; j < dimensions; ++j) {
        auto max = a.dim_max(j);
        long maxVal;
        max.foreach_piece([&maxVal](isl::set a, isl::aff expr)  {
          // this is less or equal
          maxVal = expr.get_constant_val().get_num_si() + 1;
          return isl_stat_ok;
        });
        MemRefValues[i].push_back(maxVal);
      }
      ++i;
    }
  }
  Out << "  Runner(const std::string &TraceFileName_) :\n";
  if(CreateTraceFile) {
    Out << " printer_("<< Config.RTMSize <<", TraceFileName_),";
  }
  Out << "\n";
  for(decltype(MemRefNames.size()) i = 0; i < MemRefNames.size(); ++i) {
    Out << "    " << MemRefNames[i] << "(std::make_tuple(";
    for(decltype(MemRefValues[i].size()) j = 0; j < MemRefValues[i].size(); ++j) {
      Out << MemRefValues[i][j];
      if(j < MemRefValues[i].size() - 1) {
        Out << ", ";
      }
    }
    Out << ")";
    if(CreateTraceFile) {
      // next argument is the memOffset
      Out << ", ";
      if(i == 0) {
        Out << "0";
      } else {
        Out << MemRefNames[i - 1] << ".addressAfter()";
      }
      Out << ", &printer_";
    }
    Out << ")";
    if(i < MemRefNames.size() - 1) {
      Out << ",";
    }
    Out << "\n";
  }
  Out << "  {}\n";

  std::regex DoubleBracket("\\]\\[");
  std::regex LeadingBracket("\\[");
  std::regex ClosingBracket("\\]");
  for(auto &Stmt: Program.Statements) {
    auto AstBuilder = isl::ast_build::from_context(Stmt.Domain);
    int dim = Stmt.Domain.dim(isl::dim::set);
    Out << "  void " << Stmt.Domain.get_tuple_name() << "(";
    std::string prefix = "int c";
    for(int i = 0; i < dim; ++i) {
      Out << prefix << i;
      if(i < dim - 1) {
        Out << ", ";
      }
    }
    Out << ") {\n";
    for(auto &Access: Stmt.MemoryAccesses) {
      if(Access.first.dim(isl::dim::out) == 0) {
        continue;
      }
      // LLVM_DEBUG(dbgs() << "Access: " << Access.to_str() << "\n");
      auto AccessExpr = AstBuilder.access_from(isl::pw_multi_aff::from_map(Access.first));
      auto AccessExprStr = astExprToStr(AccessExpr);
      AccessExprStr = std::regex_replace(AccessExprStr, DoubleBracket, ", ");
      if(CreateTraceFile) {
        if(Access.second) {
          AccessExprStr = std::regex_replace(AccessExprStr, LeadingBracket, ".accessAndPrint(AccessType::READ, ");
        } else {
          AccessExprStr = std::regex_replace(AccessExprStr, LeadingBracket, ".accessAndPrint(AccessType::WRITE, ");
        }

      } else {
        AccessExprStr = std::regex_replace(AccessExprStr, LeadingBracket, ".access(");
      }
      AccessExprStr = std::regex_replace(AccessExprStr, ClosingBracket, ")");
      Out << "    " << AccessExprStr << ";\n";
    }
    Out << "  }\n";

  }
  // output actual schedule
  auto AstBuilder = isl::ast_build::from_context(Program.Context);
  Out << "  void run() {\n";
  auto Node = AstBuilder.node_from_schedule_map(Schedule.intersect_domain(Program.Domains));
  Out << astNodeToStr(Node);
  Out << "  }\n";
  Out << "};\n";
  Out << "int main(int argc, char **argv) {\n";
  Out << "  std::string traceFile = \"trace\";\n"
      << "  if(argc > 1) {\n"
      << "    traceFile = std::string(argv[1]);\n"
      << "  }\n";
  Out << "  Runner runner(traceFile);\n"
      << "  runner.run();\n"
      << "  std::size_t result = 0;\n";
  // print result
  for(auto &MemRefName: MemRefNames) {
    Out << "  result += runner." << MemRefName << ".shifts;\n";
  }

  Out << "  std::cout << result << \"\\n\";\n"
    << "}\n";
}

std::string execCommand(const char* cmd);

/**
 * This calculates the number of shifts for a scop and a given schedule.
 * Currently, this is implemented by creating some C++ code, compiling it with
 * g++, and then running this, which will return the number of shifts.
 * The
 * @param S
 * @param Schedule
 * @return Precise number of iterations.
 */
std::size_t calculateNumberOfShifts(polly::Scop &S, isl::union_map Schedule, const ShiftConfig &Config = ShiftConfig());


} // hmrtm
} // polly

#undef DEBUG_TYPE
#endif //POLLY_HMSCHEDULESHIFTCALCULATOR_H
