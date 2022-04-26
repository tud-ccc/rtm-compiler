//
// Created by hauke on 06.01.20.
//

#ifndef POLLY_HMSCHEDULETREETIKZPRINTER_H
#define POLLY_HMSCHEDULETREETIKZPRINTER_H

#include "isl/isl-noexceptions.h"

#include "llvm/ADT/StringMap.h"
#include "llvm/Support/Debug.h"

#define DEBUG_TYPE "t"


static isl::set removeAllButOneDim(isl::set Set, int Dim) {
  auto NDims = Set.n_dim();
  Set = Set.remove_dims(isl::dim::set, 0, Dim);
  Set = Set.remove_dims(isl::dim::set, 1, NDims - Dim - 1);
  return Set;
}

static isl::union_map assignNewStatementNames(isl::union_set Domain) {
  std::vector<std::string> PotentialNames{"S", "T", "U", "V", "W"};
  // Only rename small kernels
  auto Renaming = isl::union_map::empty(Domain.get_space());
  auto DomainList = Domain.get_set_list();
  for(int i = 0; i < DomainList.n_set(); ++i) {
    auto RenamingPart = isl::map::identity(DomainList.get_at(i).identity().get_space());
    if(PotentialNames.size() >= DomainList.n_set()) {
      RenamingPart = RenamingPart.set_tuple_name(isl::dim::out, PotentialNames[i]);
    }
    Renaming = Renaming.unite(RenamingPart);
  }
  return Renaming;
}

static isl::union_set assignNewTupleNames(isl::union_set Domain) {
  // six names should be enough for most cases
  std::vector<std::string> PotentialNames{"i", "j", "k", "x", "y", "z"};
  auto SetList = Domain.get_set_list();
  std::vector<isl::set> Domains(SetList.n_set());
  for(int i = 0; i < SetList.n_set(); ++i) {
    Domains[i] = SetList.get_set(i);
  }
  std::sort(Domains.begin(), Domains.end(), [](const isl::set &First, const isl::set &Second) {
    return First.get_tuple_name() < Second.get_tuple_name();
  });
  llvm::StringMap<isl::set> NameRanges;
  llvm::StringMap<std::vector<std::string>> UsedNames;
  for(auto &Part: Domains) {
    auto NDims = Part.n_dim();
    auto StmtName = Part.get_tuple_name();
    UsedNames[StmtName] = std::vector<std::string>(Part.n_dim());
    auto &CurrentUsedNames = UsedNames[StmtName];
    for(unsigned i = 0; i < NDims; ++i) {
      std::string Name;
      auto Reduced = removeAllButOneDim(Part, i);
      for(auto &PotentialName: PotentialNames) {
        if(NameRanges.find(PotentialName) != NameRanges.end()
        && std::find(CurrentUsedNames.begin(), CurrentUsedNames.end(), PotentialName) == CurrentUsedNames.end()
        && Reduced.is_equal(NameRanges[PotentialName])) {
          Name = PotentialName;
          break;
        }
      }
      UsedNames[StmtName][i] = Name;
    }
    for(unsigned i = 0; i < NDims; ++i) {
      std::string Name = UsedNames[StmtName][i];
      if(Name.empty()) {
        auto j = 0;
        while(std::find(CurrentUsedNames.begin(), CurrentUsedNames.end(), PotentialNames[j]) != CurrentUsedNames.end()) {
          ++j;
        }
        Name = PotentialNames[j];
        NameRanges[Name] = removeAllButOneDim(Part, i);
        CurrentUsedNames[i] = Name;
      }
    }
  }
  Domain = isl::union_set::empty(Domain.get_space());
  for(auto &Part: Domains) {
    auto NDims = Part.n_dim();
    auto StmtName = Part.get_tuple_name();
    auto NameList = UsedNames[StmtName];
    for(unsigned i = 0; i < NDims; ++i) {
      isl_id * LoopIdRaw = isl_id_alloc(Part.get_ctx().get(), NameList[i].c_str(), nullptr);
      auto LoopId = isl::manage(LoopIdRaw);
      Part = Part.set_dim_id(isl::dim::set, i, LoopId);
    }
    Domain = Domain.unite(Part);
  }
  return Domain;
}

static std::string escapeForLatex(const std::string &Str) {
  std::regex OpenBrace("\\{");
  std::regex ClosingBrace("\\}");
  auto Res = std::regex_replace(Str, OpenBrace, "\\{");
  Res = std::regex_replace(Str, ClosingBrace, "\\}");
  return Res;
}

template<typename Stream>
static void printWhitespace(Stream &Out, int Whitespace) {
  for(int i = 0; i < Whitespace; ++i) {
    Out << " ";
  }
}

template<typename Stream>
static void printScheduleNodeToStream(const isl::schedule_node &Node,
    Stream &Out, int Whitespace, bool PrintLeaves) {
  auto NodeType = isl_schedule_node_get_type(Node.get());
// isl_schedule_node_band,
//	isl_schedule_node_context,
//	isl_schedule_node_domain,
//	isl_schedule_node_expansion,
//	isl_schedule_node_extension,
//	isl_schedule_node_filter,
//	isl_schedule_node_leaf,
//	isl_schedule_node_guard,
//	isl_schedule_node_mark,
//	isl_schedule_node_sequence,
//	isl_schedule_node_set
  Out << "node { ";
  switch (NodeType) {
    case isl_schedule_node_band: {
      Out << "band: $";
      auto Sched = isl::manage(isl_schedule_node_band_get_partial_schedule_union_map(Node.get()));
      Out << escapeForLatex(Sched.to_str());
      Out << "$";
      break;
    }
    case isl_schedule_node_domain: {
      Out << "domain: $";
      auto Domain = isl::manage(isl_schedule_node_domain_get_domain(Node.get()));
      Out << escapeForLatex(Domain.to_str());
      Out << "$";
      break;
    }
    case isl_schedule_node_filter: {
      Out << "filter: $";
      auto Filter = isl::manage(isl_schedule_node_filter_get_filter(Node.get()));
      Out << escapeForLatex(Filter.to_str());
      Out << "$";
      break;
    }
    case isl_schedule_node_sequence: {
      Out << "sequence";
      break;
    }

    case isl_schedule_node_leaf: {
      Out << "leaf";
      break;
    }
    default:
      Out << "ignored";
      break;
  }
  Out << " }\n";
  for(int i = 0; i < Node.n_children(); ++i) {
    auto Child = Node.get_child(i);
    if(isl_schedule_node_get_type(Child.get()) != isl_schedule_node_leaf
    || PrintLeaves) {
      printWhitespace(Out, Whitespace);
      Out << "child { ";
      printScheduleNodeToStream(Child, Out, Whitespace + 2, PrintLeaves);
      printWhitespace(Out, Whitespace);
      Out << "}\n";
    }
  }
}

template<typename Stream>
static void printScheduleToStream(isl::schedule Schedule, Stream &Out) {
  auto Mapping = assignNewStatementNames(Schedule.get_domain());

  Out << Mapping.to_str() << "\n";
  Out << assignNewTupleNames(Schedule.get_domain()).to_str() << "\n";
  auto Root = Schedule.get_root();
  Out << "\\";
  printScheduleNodeToStream(Root, Out, 2, true);
  Out << ";";
}
#undef DEBUG_TYPE
#endif //POLLY_HMSCHEDULETREETIKZPRINTER_H
