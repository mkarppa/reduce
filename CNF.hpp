#ifndef REDUCE2_CNF_HPP
#define REDUCE2_CNF_HPP

#include <vector>
#include <cstdlib>
#include <iostream>
#include <set>
#include <map>

namespace reduce {
  class CNF {
  public:
    // variable type (technically unsigned)
    typedef int Variable;
    // literal type (signed)
    typedef int Literal;
    // clause indexing
    typedef int ClauseIdx;
    // clause type (vector of literals)
    typedef std::vector<Literal> Clause;

    int getVariableCount() const {
      return variables.size();
    }

    int getClauseCount() const {
      return clauses.size();
    }

    void addClause(const Clause& clause) {
      for (Variable v : clause)
        if (!variables.count(abs(v)))
          throw std::runtime_error("Attempted to add an invalid clause!");
      clauses.push_back(clause);
    }

    CNF() = default;
    
    /**
     * Initialize CNF from a DIMACS stream
     */
    explicit CNF(std::istream& is) {
      parseDimacs(is);
    }

    bool empty() const {
      return variables.size() == 0 && clauses.size() == 0;
    }

    std::vector<Clause>::const_iterator begin() const {
      return clauses.begin();
    }

    std::vector<Clause>::const_iterator end() const {
      return clauses.end();
    }

    const Clause& operator[](ClauseIdx i) const {
      return clauses[i];
    }

    const std::map<Literal,std::vector<ClauseIdx> >& getLiteralOccurrences()
      const {
      return literalOccurrences;
    }


  private:
    void parseDimacs(std::istream& is);
    
    std::set<Variable> variables;
    std::vector<Clause> clauses;
    // which clauses contain which literals
    std::map<Literal,std::vector<ClauseIdx> > literalOccurrences;
  };
}

#endif // REDUCE2_CNF_HPP
