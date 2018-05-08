#ifndef REDUCE2_REDUCER_HPP
#define REDUCE2_REDUCER_HPP

#include "ConflictChecker.hpp"
#include "Stack.hpp"
#include "VariableMapping.hpp"
#include "CNF.hpp"
#include "common.hpp"
#include <vector>
#include <deque>
#include <memory>

namespace reduce {
  /*************************************************** Reducer data structure. */

  class Reducer {
  public:
    // assumes responsibility for deleting workStack
    Reducer(bool verbose, bool printAutSizes, bool checkConflicts,
            const reduce::CNF& inCnf,
            const std::vector<int>& prefix, 
            long threshold, const Graph& graph, const VariableMapping& varMap,
            Stack* workStack, const std::vector<double>& estimationProbabilities = std::vector<double>()
#ifdef WITH_OPEN_MPI
            , const std::vector<std::vector<int> >& lvlsForSlaves = std::vector<std::vector<int> >()
#endif // WITH_OPEN_MPI
            );
    

    
    // disable copy and move operations
    Reducer(const Reducer&) = delete;
    Reducer& operator=(const Reducer&) = delete;
    Reducer(Reducer&&) = delete;
    Reducer& operator=(Reducer&&) = delete;

    // typedef const int* Assignment;
    struct Assignment {
      std::vector<int> vars;
      std::vector<int> vals;
      int autSize;
    };
    
    /**
     * Print a prefix assignment obtained from the reducer.
     */
    void printAssignment(std::ostream& out, const Assignment& a) const;


    void printCnf(std::ostream& out, 
                  const std::string& fmt, 
                  int header_var_adjust, 
                  int header_clause_adjust) const;

    const CNF& getClauses() const {
      return clauses;
    }

    const Graph& getBaseGraph() const {
      return base;
    }

    const VariableMapping& getVarMap() const {
      return varMap;
    }

    size_t getPrefixLength() const {
      return prefixSequence.size();
    }

    void printStats(std::ostream& out, int level) const;

    void getPrefixAssignments(std::vector<Assignment>& out);

    inline size_t getConflictCount() const {
      return conflictCount;
    }
    
  private:
    CNF clauses; // The clauses 
    Graph base; // The base graph
    VariableMapping varMap;
    ConflictChecker conflictChecker;
    bool verbose;         /* Verbose output? */
    bool print_aut_sizes; /* Print automorphism group sizes wrt cubes? */
    bool checkConflicts; // Whether to check for conflicts in CNF case
    long autSizeThreshold;               /* Threshold size for automorphism group. */
    std::vector<int> prefixSequence; // The prefix sequence
    std::vector<size_t> statGen; /* Generated assignments. */
    std::vector<size_t> statCan; /* Canonical assignments. */
    std::vector<size_t> statOut; /* Assignments output. */
    std::vector<std::vector<int> > orbits; /* Indicators for prefix element orbits.*/
    std::vector<std::vector<int> > travInd; /* Traversal indicators. */
    std::vector<std::vector<Permutation> > traversals; /* The traversal permutations. */  
    
    std::unique_ptr<Stack> work;
    std::vector<int> scratch; // scratch

    std::vector<double> estimationProbabilities;

    size_t conflictCount;

#ifdef WITH_OPEN_MPI
    std::vector<int> mpiBuffer; // buffer for message passing
    int mpiRank;
    int mpiSize;
    std::vector<std::vector<int> > lvlsForSlaves;
#endif // WITH_OPEN_MPI
    
    void initialize();
    void buildPrefix();
    bool stacksEqual() const;
    const int* getPrefixAssignment();
    void initWork();

    // push next variable value pairs
    void pushNext(int* seedMin, int* vars, int* vals, int size);
    
    void getPrefixAssignmentsSolo(std::vector<Assignment>& out);
#ifdef WITH_OPEN_MPI
    void getPrefixAssignmentsMaster(std::vector<Assignment>& out);
    void getPrefixAssignmentsSlave();
#endif // WITH_OPEN_MPI
  };
}

#endif // REDUCE2_REDUCER_HPP
