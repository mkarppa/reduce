#ifndef REDUCE2_CONFLICTCHECKER_HPP
#define REDUCE2_CONFLICTCHECKER_HPP

#include "VariableMapping.hpp"
#include "CNF.hpp"

namespace reduce {
  // CNF conflict checker
  class ConflictChecker {
  public:
    explicit ConflictChecker(const CNF& cnf,
                             const VariableMapping* varMap = nullptr);

    // returns true if the given partial assignment of k variables causes a
    // conflict (ie. is unsatisfiable as such)
    bool hasConflict(const int* vars, const int* vals, size_t k) const;
    
  private:
    std::map<CNF::Literal,std::vector<CNF::ClauseIdx> > literalOccurrences;
    std::vector<size_t> clauseSizes;
    const VariableMapping* varMap;
  };
}

#endif // REDUCE2_CONFLICTCHECKER_HPP
