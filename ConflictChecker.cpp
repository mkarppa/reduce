#include "ConflictChecker.hpp"
#include <cassert>

using std::vector;
using std::cerr;
using std::endl;

namespace reduce {
  ConflictChecker::ConflictChecker(const CNF& cnf,
                                   const VariableMapping* varMap) :
    literalOccurrences(cnf.getLiteralOccurrences()),
    clauseSizes(cnf.getClauseCount()), varMap(varMap) {
    auto it = clauseSizes.begin();
    for (auto jt = cnf.begin(); jt != cnf.end(); ++jt)
      *it++ = jt->size();
  }

  
  
  bool ConflictChecker::hasConflict(const int* vars, const int* vals, size_t k) const {
    vector<size_t> satisfiableVariables(clauseSizes);
    for (size_t i = 0; i < k; ++i) {
      int var = varMap ? varMap->graphToCnf(vars[i])+1 : vars[i];
      int val = vals[i];
      assert(val == 1 || val == 0);
      CNF::Literal l = (1-2*val)*var; // take the anti-literal
      assert(literalOccurrences.count(l));
      for (auto c : literalOccurrences.find(l)->second) {
        if (--satisfiableVariables[c] == 0)
          return true;
      }
    }
    return false;
  }
}
