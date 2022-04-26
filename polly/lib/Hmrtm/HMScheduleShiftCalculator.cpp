//
// Created by hauke on 16.10.19.
//

#include "polly/HMScheduleShiftCalculator.h"

#define DEBUG_TYPE "polly-hm-trace-generation"

using namespace polly::hmrtm;

polly::hmrtm::Statement::Statement(const isl::set &Domain)
    : Domain(Domain), MemoryAccesses() {}
polly::hmrtm::Statement::Statement(const Statement &s)
    : Domain(s.Domain), MemoryAccesses(s.MemoryAccesses) {
  //LLVM_DEBUG(dbgs() << "Copied statement!\n");
}
void polly::hmrtm::Statement::addAccess(const isl::map &Access, bool IsRead) {
  MemoryAccesses.push_back(std::make_pair(Access, IsRead));
}
polly::hmrtm::ShiftCalculationProgram::ShiftCalculationProgram(const isl::union_set &Domains, const isl::set &Context)
    : Domains(Domains), Context(Context), Statements(), Variables() {}
void polly::hmrtm::ShiftCalculationProgram::addStatement(polly::hmrtm::Statement &&Stmt) {
  for(auto &Access: Stmt.MemoryAccesses) {
    Variables.insert(Access.first.get_tuple_name(isl::dim::out));
  }
  Statements.push_back(std::move(Stmt));
}
polly::hmrtm::ShiftCalculationProgram polly::hmrtm::createShiftProgram(polly::Scop &S) {
  // domains need to match
  //assert(S.getDomains().is_equal(schedule.domain()));
  //S.getDomains().dump();
  //Schedule.domain().dump();
  auto LoopSize = 96;
  auto Context = S.getContext();
  auto LSp = isl::local_space(Context.get_space());
  auto Params = Context.dim(isl::dim::param);
  // first assign values to param
  for(unsigned i = 0; i < Params; ++i) {
    auto Constraint = isl::constraint::alloc_equality(LSp);
    Constraint = Constraint.set_constant_si(LoopSize);
    Constraint = Constraint.set_coefficient_si(isl::dim::param, i, -1);
    Context = Context.add_constraint(Constraint);
  }
  LLVM_DEBUG(dbgs() << "Context: " << Context.to_str());

  //auto AstBuilder = isl::ast_build::from_context(Context);
  // insert parameter adjustment here
  auto Domains = S.getDomains().gist_params(Context);
//  Domains.dump();
//  std::string NewParamBaseName = "POLYHEDRAL_SIZE_PARAM_";
//  auto GenerateNextName = [&NewParamBaseName]() {
//    NewParamBaseName = NewParamBaseName + std::string("I");
//    return NewParamBaseName;
//  };
//  llvm::StringMap<std::string> paramReplacement;
  ShiftCalculationProgram SCP(Domains, Context);
  for (ScopStmt &StmtScop : S) {
    Statement Stmt(StmtScop.getDomain().gist_params(Context));
    auto AstBuilder = isl::ast_build::from_context(Stmt.Domain);
    for (MemoryAccess *MA : StmtScop) {
      auto MemAccessMap = MA->getLatestAccessRelation();
      //AstBuilder.access_from(isl::pw_multi_aff::from_map(MemAccessMap)).dump();
      Stmt.addAccess(MemAccessMap, MA->getType() == MemoryAccess::AccessType::READ);
    }
    SCP.addStatement(std::move(Stmt));
  }
  return SCP;
}
std::string polly::hmrtm::astExprToStr(const isl::ast_expr &ast) {
  auto printer = isl_printer_to_str(ast.get_ctx().get());
  printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
  printer = isl_printer_print_ast_expr(printer, ast.get());
  std::string out(isl_printer_get_str(printer));
  isl_printer_free(printer);
  return out;
}
std::string polly::hmrtm::astNodeToStr(const isl::ast_node &ast) {
  auto printer = isl_printer_to_str(ast.get_ctx().get());
  printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
  printer = isl_printer_print_ast_node(printer, ast.get());
  std::string out(isl_printer_get_str(printer));
  isl_printer_free(printer);
  return out;
}
std::string polly::hmrtm::execCommand(const char *cmd) {
  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
    llvm_unreachable("Temp program must exist!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  return result;
}
std::size_t polly::hmrtm::calculateNumberOfShifts(polly::Scop &S,
                                                  isl::union_map Schedule,
                                                  const ShiftConfig &Config) {

  // Cannot calculate the shifts in this case!
  if(Schedule.dim(isl::dim::param) != 0) {
    return 0;
  }
  LLVM_DEBUG(Config.print(dbgs()));
  std::string OutputFilePrefix = (Config.OutputDirectory + "/trace__"
      + Config.Position + "__"+ Config.KernelName + "__"
      + S.getFunction().getName() + "__" + S.getNameStr()).str();
  auto SCP = polly::hmrtm::createShiftProgram(S);
  SmallVector<char, 30> path;
  SmallVector<char, 30> execPath;
  int fd;
  llvm::sys::fs::createTemporaryFile("Hallo", "cpp", fd, path);
  llvm::sys::fs::createUniquePath("shifts-calculator-%%%%%%%%%%", execPath, true);
  llvm::raw_fd_ostream tempSourceFile(fd, true);
  //printProgramToStream(dbgs(), SCP);
  printProgramToStream(tempSourceFile, SCP, Schedule, Config, OutputFilePrefix + ".nvt");
  tempSourceFile.flush();
  std::string arg = ("c++ -O3 -DNDEBUG " + path + " -o " + execPath).str();
  auto compilationResult = std::system(arg.c_str());
  if(compilationResult != 0) {
    llvm_unreachable("Call to c++ failed!");
  }
  arg = (Twine(execPath) + " \"" + OutputFilePrefix + ".nvt\"").str();
  LLVM_DEBUG(dbgs() << Twine(path).str() << ": " << arg);
  auto result = execCommand(arg.c_str());
  std::stringstream sstream(result);
  std::size_t shifts;
  sstream >> shifts;
  LLVM_DEBUG(dbgs() << "Got shifts: " << shifts << "\n");
  tempSourceFile.close();
  if(!Config.OutputDirectory.empty() && !Config.KernelName.empty()) {
    LLVM_DEBUG(dbgs() << "Copying files...\n");
    llvm::sys::fs::copy_file(path, OutputFilePrefix  + ".cpp");
    llvm::sys::fs::copy_file(execPath, OutputFilePrefix + ".exec");

  }
  llvm::sys::fs::remove(path);
  llvm::sys::fs::remove(execPath);

  return shifts;
}

#undef DEBUG_TYPE