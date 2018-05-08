#include "Reducer.hpp"
#include "mpiaux.hpp"
#include "gmp.h"
#include <random>
#include <cstring>
#include <cassert>

using std::to_string;
using std::runtime_error;
using std::vector;
using std::cerr;
using std::endl;
using std::string;

namespace reduce {

  /***************************************** Subroutines for orbit traversals. */

  static vector<Permutation> traversal_prepare(int root, const Graph& g) {
    int n = g.getOrder();
    if (root < 0 || root >= n)
      throw runtime_error("bad root");
    vector<int> ind(n,-1);
    vector<Graph::Vertex> list = g.orbitOf(root);
    int rootpos = std::find(list.begin(), list.end(), root) - list.begin();
    for (size_t i = 0; i < list.size(); ++i)
      ind[list[i]] = i;
       
    vector<Permutation> t(list.size());
    ind[root] = -1;
    t[rootpos] = Permutation(n);

    // not sure if this is very stylish
    // anyway, repeat while there is a non-negative indicator left
    while([&]() {
        for (auto l : list) {
          if(ind[l] >= 0) {
            return true;
          }
        }
        return false;
      }()) {      
      for (auto p : g.getAutGenPerm()) {
        for (size_t j = 0; j < list.size(); j++) {
          int u = list[j]; // choose a vertex
          int v = p(u); // permute it
          int q = ind[v];  // check if the resulting vertex has already been visited
          if(q >= 0 && ind[u] < 0) {
            // if not, record the permutation
            t[q] = p * t[j];
            ind[v] = -1; // mark the vertex as visited
          }
        }
      }
    }

    for (size_t j = 0; j < list.size(); ++j)
      if (t[j](root) != static_cast<int>(list[j]))
        throw runtime_error("bad traversal");

    return t;
  }

  
  /********************************* Get a prefix assignment from the reducer. */

  static int autOrderTrunc(const Graph& g) {
    mpz_class autOrder = g.getAutOrder();
    int autTrunc = 999999999;
    if (cmp(autOrder,autTrunc) < 0)
      autTrunc = autOrder.get_si();
    return autTrunc;
  }

  /****************************************** Initialize a configured reducer. */

  // compute indicators for minimum elements on orbits
  static vector<int> orbitMinimumIndicators(const Graph& g,
                                            int *relabel) {
    int n = g.getOrder();
    const vector<int>& p = g.orbitCells();
    const vector<int>& c = g.getOrbits();
    vector<int> ind(n,0);
    for (int i = 0; i < n; i++) {
      if (relabel != nullptr)
        ind[relabel[p[i]]] = 1;
      else
        ind[p[i]] = 1;
      int j = i+1;
      while (j < n && c[p[i]] == c[p[j]])
        ++j;
      i = j-1;
    }
    return ind;
  }


  
  /*
    print the automorphism group order with a given prefix and postfix string
  */
  static void prettyPrintAutOrder(std::ostream& out,
                                  const Graph& g,
                                  const string& prefix,
                                  const string& postfix) {
    out << prefix << g.getAutOrder().get_str() << postfix;
  }


  
  static void printAutOrder(std::ostream& out, const Graph& g) {
    prettyPrintAutOrder(out,g,"\n   |Aut| = ","");
  }


  
  static void print_orbit_perms(FILE *out, const Graph& g, const VariableMapping& varMap)
  {
    // this function is ugly but not important enough for reimplementation
    int n = g.getOrder();
    int l = varMap.getVariableCount();
    const vector<int>& p = g.orbitCells();
    const vector<int>& c = g.getOrbits();
    vector<int> q(n,0);
    for(int u = 0; u < l; u++)
      q[varMap.getVariable(u)] = 1;

    for(int i = 0; i < n; i++) {
      int j = i+1;
      for(; j < n && c[p[i]] == c[p[j]]; j++)
        ;
      if(q[p[i]] > 0) {
        fprintf(out, "orbit: ");
        print_int_array(out, j - i, &p[0] + i);
        fprintf(out, "\n");
        for (auto a : g.getAutGen()) {
          fprintf(out, "       ");
          for(int u = 0; u < j - i; u++)
            q[p[u + i]] = 2;
          int u = 0;
          int num_fixed = 0;
          int num_moved = 0;
          while(u < l) {
            int z = varMap.getVariable(u);
            int first = 1;
            int len = 0;
            if(q[z] == 2) {
              int w = z;
              do {                        
                q[w] = 1;
                fprintf(out, 
                        "%s%d",
                        first == 1 ? "(" : " ",
                        w + 1);
                first = 0;
                w = a[w];
                fprintf(out, "%s",
                        w == z ? ")" : "");
                len++;
              } while(z != w);
              if(len == 1)
                num_fixed += len;
              if(len >= 2)
                num_moved += len ;                     
            }
            u++;
          }
          fprintf(out, " -- fix = %d, move = %d\n", num_fixed, num_moved);
        }
      }
      i = j-1;
    }
  }



  void Reducer::printStats(std::ostream& out, int l) const {
    char temp[64];
    sprintf(temp, 
            "c %7d %14ld %14ld %14ld", 
            l+1, 
            statGen[l], 
            statCan[l],
            statOut[l]);
    out << temp << endl;
  }


  
  static void printGraphStats(const Graph& g, int idx, const VariableMapping& varMap) {
    cerr << "graph [" << idx << "]:";
    printAutOrder(cerr, g);
    cerr << endl;
    cerr << "   orbits = [";
    g.printOrbits(cerr, varMap);
    cerr << "]" << endl;
  }


  
  void Reducer::buildPrefix() {
    assert(orbits.size() == prefixSequence.size());
    
    Graph g(base);
    for(size_t k = 0; k < prefixSequence.size(); ++k) {
      push_time();

      assert(orbits[k].size() == static_cast<size_t>(base.getOrder()));

      cerr << "graph [" << k << "]:";
      printAutOrder(cerr, g);
      cerr << endl;

      if (k == 0) {
        /* Check the base graph against the variable and value lists. */
        const vector<int>& p = g.orbitCells();
        const vector<int>& c = g.getOrbits();

        vector<int> q(g.getOrder(), 0);
        for (size_t j = 0; j < varMap.getVariableCount(); ++j)
          q[varMap.getVariable(j)] = 1; // flag variables
      
        for (Graph::Size s = 0; s < g.getOrder(); ++s) {
          Graph::Size u = s+1;
        
          // while the span belongs to the same orbit
          while (u < g.getOrder() && c[p[s]] == c[p[u]])
            ++u;
        
          for (Graph::Size j = s+1; j < u; ++j) {
            if(q[p[j]] != q[p[s]]) {
              throw runtime_error("variable list is not a union of "
                                  "orbits of base graph "
                                  "(" + to_string(p[j] + 1) +
                                  " and " + to_string(p[s]+1) +
                                  " have different orbits)");
            }
          }
          s = u-1;
        }

        for (Graph::Size j = 0; j < g.getOrder(); ++j)
          for (Graph::Size s = 0;
               static_cast<size_t>(s) < varMap.getValueCount(); ++s)
            if (p[j] == static_cast<int>(varMap.getValue(s)) &&
               ((j > 0 && c[p[j-1]] == c[p[j]]) ||
                (static_cast<int>(j) < static_cast<int>(g.getOrder()) - 1 &&
                                     c[p[j]] == c[p[j+1]])))
              throw runtime_error("value vertex (" +
                                  to_string(varMap.getValue(s) + 1) +
                                  ") is not fixed by the automorphism "
                                  "group of the base graph");
      }
      cerr << "   orbits = [";
      g.printOrbits(std::cerr, varMap);
      cerr << "]" << endl;

      if (verbose) 
        print_orbit_perms(stderr, g, varMap);
    
      cerr << "prefix[" << (k+1) << "] = " << (prefixSequence[k]+1) << ":";

      push_time();
      traversals[k] = traversal_prepare(prefixSequence[k], g);
      pop_print_time("traversal");

      vector<int> a(traversals[k].size());
      for (size_t j = 0; j < traversals[k].size(); ++j)
        a[j] = traversals[k][j](prefixSequence[k]);
      for (Graph::Size i = 0; i < base.getOrder(); ++i)
        travInd[k][i] = 0;
      for (size_t j = 0; j < traversals[k].size(); ++j)
        travInd[k][a[j]] = 1;
      fprintf(stderr, "\n   traversal: ");
      print_int_array(stderr, traversals[k].size(), &a[0]);
      fprintf(stderr, " [length = %lu]\n", traversals[k].size());

      g.addEdge(prefixSequence[k], varMap.getValue(0));
      for (Graph::Size j = 0; j < base.getOrder(); ++j)
        orbits[k][j] = g.sameOrbit(prefixSequence[k], j);

      pop_print_time("prefix_total");
      fprintf(stderr, "\n");
    }

    if (prefixSequence.size() > 0) {
      printGraphStats(g,prefixSequence.size(),varMap);
      if(verbose) 
        print_orbit_perms(stderr, g, varMap);
    }
  }

  

  void Reducer::initialize() {
    /* Test prefix for repeated elements. */
    vector<int> q(prefixSequence);
    std::sort(q.begin(), q.end());
    for (size_t i = 1; i < prefixSequence.size(); ++i) {
      if(q[i-1] == q[i])
        throw runtime_error("prefix repeats an element (" +
                            to_string(q[i] + 1) + ")");
    }

    for (int i : prefixSequence) {
      if (varMap.graphToCnf(i) == static_cast<int>(VariableMapping::NOT_A_VARIABLE))
        throw runtime_error("prefix element (" + to_string(i+1) +
                            ") is not a declared variable vertex");
    }


    push_time();   

    buildPrefix();
   
    cerr << "init:";
    pop_print_time("reducer_initialize");
    cerr << endl;
  }



  // generate a value q uniformly at random in [0,1] and return true if q < p
  static bool coin(double p) {
    static std::mt19937 rng(std::chrono::system_clock::now().time_since_epoch().count());
    static std::uniform_real_distribution<double> distribution(0.,1.);
    return distribution(rng) < p;
  }
  

  const int* Reducer::getPrefixAssignment() {
    int n = base.getOrder();
    int d = varMap.getValueCount();

    if (prefixSequence.size() == 0) // nothing to do
      return nullptr;

    // work initialization should be done before calling this function   

    StackElement topElem;
    while (work->pop(topElem)) {
      /* Pop the stack top. */

      int size = topElem.getSize(); 
      int* vars = topElem.getVars();
      int* vals = topElem.getVals();
     
      /* Find the current variable. */
      int lvl = size - 1;
      int current = -1;
      int current_idx = -1;
      for (size_t j = 0; j < traversals[lvl].size(); ++j) {
        for (int i = 0; i < size; i++) {             
          if (vars[i] == traversals[lvl][j](prefixSequence[lvl])) {
            current = j;
            current_idx = i;
          }
        }
      }
      if (current == -1)
        throw runtime_error("no current variable");

      int current_val = vals[current_idx];
     
      if (current_val < d) {
        ++statGen[lvl];
        
        // if we are dealing with CNF and have enabled conflict
        // checking, check for conflicts and move up in the stack if a
        // conflict is found
        if (!clauses.empty() && checkConflicts) {
          if (conflictChecker.hasConflict(vars, vals, size)) {
            if (verbose) {
              cerr << "Conflict found: ";
              for (int i = 0; i < size; ++i)
                cerr << (varMap.graphToCnf(vars[i])+1)*(2*vals[i]-1)
                     << (i == size-1 ? "\n" : " ");
            }
            ++conflictCount;
            continue; // move to the next element in case of conflict
          }
        }



        /* Process stack top. */
        vector<int> nu(n);
        for(int i = 0; i < n; i++)
          nu[traversals[lvl][current](i)] = i;
        if(nu[vars[current_idx]] != prefixSequence[lvl])
          throw runtime_error("bad nu");
        Graph g(base);
        for (int i = 0; i < size; i++) {
          if (i != current_idx) {
            g.addEdge(vars[i], varMap.getValue(vals[i]));
          }
          else {
            g.addEdge(vars[i], varMap.getValue(current_val));
          }
        }
        const vector<int>& lab = g.canLab();
        int qlab = -1;
        int t = 0;          
        for (; t < n; t++) {
          qlab = lab[t];
          if (orbits[lvl][nu[qlab]])
            break;
        }
        if (t == n)
          throw runtime_error("bad qlab");
        
        if (g.sameOrbit(qlab, vars[current_idx])) {
          /* Top was accepted by isomorph rejection. */
          ++statCan[lvl];
                
          /* Normalize top to scratch. */
          scratch[0] = size;
          int *norm_vars = &scratch[1]; 
          int *norm_vals = &scratch[size+1]; 
          for (int i = 0; i < size; i++) {
            norm_vars[i] = nu[vars[i]];
            if (i != current_idx) {
              norm_vals[i] = vals[i];
            }
            else {
              norm_vals[i] = current_val;
            }
          }

          int aut = autOrderTrunc(g);
          scratch[2*size+1] = aut;

          if (size == static_cast<int>(prefixSequence.size()) ||
              aut <= autSizeThreshold) {
            for (int i = 0; i < size; i++) 
              norm_vals[i] = varMap.getValue(norm_vals[i]);
            if (print_aut_sizes) 
              prettyPrintAutOrder(cerr, g, "a   |Aut| = ", "\n");

            /* Report to caller. */
            ++statOut[lvl];
            return &scratch[0];
          }
          else {
            vector<int> expVars(size+1);
            vector<int> expVals(size+1);
                    
            for(int i = 0; i < size; i++) {
              expVars[i] = norm_vars[i];
              expVals[i] = norm_vals[i];
            }

            /* Save minima of (normalised) automorphism orbits. */
            vector<int> seedMin = orbitMinimumIndicators(g, &nu[0]);

            /* First var is minimum in its seed-automorphism orbit. */
            int s = 0;
            for (; s < static_cast<int>(traversals[lvl+1].size()); s++) {
              if (seedMin[traversals[lvl+1][s](prefixSequence[lvl+1])]) {
                expVars[size] = 
                  traversals[lvl+1][s](prefixSequence[lvl+1]);
                break;
              }
            }
            if (s == static_cast<int>(traversals[lvl+1].size()))
              throw runtime_error("no minimum found in extending orbit");

            if (estimationProbabilities.size() == 0 ||
                coin(estimationProbabilities[lvl]))
              pushNext(&seedMin[0], &expVars[0], &expVals[0], size+1);            
          }
        }
      }
    }
    if (!work->empty())
      throw runtime_error("work stack out of balance");
    /* Work done. */

    return nullptr;
  }


  
  void Reducer::pushNext(int* seedMin, int* vars,
                         int* vals, int size) {
    int d = varMap.getValueCount();
    int lvl = size - 1;
    int current = -1;
    int current_idx = -1;
    for (size_t j = 0; j < traversals[lvl].size(); ++j) {
      for (int i = 0; i < size; i++) {             
        if (vars[i] == traversals[lvl][j](prefixSequence[lvl])) {
          current = j;
          current_idx = i;
        }
      }
    }
    if (current == -1)
      throw runtime_error("no current variable");
    
    vector<int> seedMinimumVars;
    for(int p = current; p < static_cast<int>(traversals[lvl].size()); ++p) {
      if(seedMin[traversals[lvl][p](prefixSequence[lvl])]) {
        seedMinimumVars.push_back(traversals[lvl][p](prefixSequence[lvl]));
      }
    }

    /* Next variable must be minimum in its seed-automorphism orbit. */
    for (auto it = seedMinimumVars.rbegin(); it != seedMinimumVars.rend(); ++it) {
      int p = *it;
      vars[current_idx] = p; 
      for (int r = d; r >= 0; --r) {
        /* Save next value */       
        vals[current_idx] = r;
        work->push(StackElement(vector<int>(vars,vars+size),
                                vector<int>(vals,vals+size),
                                vector<int>(seedMin,seedMin+base.getOrder())));
      }
    }
  }


  
  /********************** Print a prefix assignment obtained from the reducer. */
  void Reducer::printAssignment(std::ostream& out, const Reducer::Assignment& a) const {
    int size = a.vars.size();
    assert(a.vars.size() == a.vals.size());
    const int *vars = &a.vars[0]; 
    const int *vals = &a.vals[0]; 
    for(int i = 0; i < size; i++) {
      size_t j, jj;
      for(j = 0; j < varMap.getVariableCount(); j++) {
        if (vars[i] == static_cast<int>(varMap.getVariable(j)))
          break;          
      }
      for(jj = 0; jj < varMap.getValueCount(); jj++) {
        if (vals[i] == static_cast<int>(varMap.getValue(jj)))
          break;          
      }
      if (j == varMap.getVariableCount() || jj == varMap.getValueCount())
        throw std::runtime_error("no data for assignment");
      out << varMap.variableLegend(j) << " -> " << varMap.valueLegend(jj)
          << (i == size-1 ? "\n" : ", ");
    }
  }



  void Reducer::printCnf(std::ostream& out, 
                         const std::string& fmt, 
                         int header_var_adjust, 
                         int header_clause_adjust) const {
    if(clauses.empty())
      throw std::runtime_error("do not have CNF to print");

    bool do_header_counts = true;
    if(header_var_adjust < 0 ||
       header_clause_adjust < 0)
      do_header_counts = false;

    if(do_header_counts) {
      out << "p " << fmt << " "
          << clauses.getVariableCount() + header_var_adjust
          << " "
          << clauses.getClauseCount() + header_clause_adjust
          << std::endl;
    }
    else {
      out << "p " << fmt << std::endl;
    }
    
    for (const CNF::Clause& c : clauses) {
      bool first = true;
      for (CNF::Literal l : c) {
        out << (first ? "" : " ") << l;
        first = false;
      }
      out << (first ? "" : " ") << "0" << std::endl;
    }
  }



  Reducer::Reducer(bool verbose, bool printAutSizes, bool checkConflicts,
                   const reduce::CNF& inCnf,
                   const std::vector<int>& prefix, 
                   long threshold, const Graph& graph,
                   const VariableMapping& varMap,                   
                   Stack* workStack,
                   const std::vector<double>& estimationProbabilities
#ifdef WITH_OPEN_MPI
            , const std::vector<std::vector<int> >& lvlsForSlaves
#endif // WITH_OPEN_MPI
                   ) :
    clauses(inCnf),
    base(graph),
    varMap(varMap),
    conflictChecker(inCnf, &this->varMap),
    verbose(verbose),
    print_aut_sizes(printAutSizes),
    checkConflicts(checkConflicts),
    autSizeThreshold(threshold),
    prefixSequence(prefix),
    statGen(prefix.size(), 0),
    statCan(prefix.size(), 0),
    statOut(prefix.size(), 0),
    orbits(vector<vector<int> >(prefixSequence.size(),
                                vector<int>(base.getOrder()))),
    travInd(vector<vector<int> >(prefixSequence.size(),
                                 vector<int>(base.getOrder()))),
    traversals(prefixSequence.size()),
    work(workStack),
    scratch(2*prefix.size()+2),
    estimationProbabilities(estimationProbabilities),
    conflictCount(0)
#ifdef WITH_OPEN_MPI
    , mpiBuffer((2*prefix.size()+2+base.getOrder())*(prefix.size()+1)),
    mpiRank(-1),
    mpiSize(-1),
    lvlsForSlaves(lvlsForSlaves)
#endif // WITH_OPEN_MPI
  {
#ifdef WITH_OPEN_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif // WITH_OPEN_MPI

    assert(estimationProbabilities.size() == 0 ||
           estimationProbabilities.size() == prefixSequence.size());
    
    std::cerr << "input: n = " << base.getOrder()
              << ", m = " << base.getNumEdges()
              << ", v = " << varMap.getVariableCount()
              << ", r = " << varMap.getValueCount()
              << ", k = " << prefixSequence.size()
              << ", t = " << autSizeThreshold
              << std::endl;
    std::cerr << "prefix:";
    for (auto p : prefixSequence)
      cerr << " " << p;
    cerr << endl;
    pop_print_time("reducer_parse");
    std::cerr << std::endl;

    initialize();
  }



  void Reducer::getPrefixAssignments(std::vector<Assignment>& out) {
#ifdef WITH_OPEN_MPI
    if (mpiSize == 1 && mpiRank == 0)
#endif // WITH_OPEN_MPI
      getPrefixAssignmentsSolo(out);
#ifdef WITH_OPEN_MPI
    else if (mpiSize > 1 && mpiRank == 0)
      getPrefixAssignmentsMaster(out);
    else if (mpiSize > 1 && mpiRank > 0)
      getPrefixAssignmentsSlave();
    else
      throw runtime_error("Invalid MPI rank (" + to_string(mpiRank) +
                          ") or size (" + to_string(mpiSize) + ")");
#endif // WITH_OPEN_MPI
  }
  


  void Reducer::getPrefixAssignmentsSolo(vector<Reducer::Assignment>& assignments) {
    initWork();
    assignments = vector<Assignment>();
    const int* a = nullptr;
    while ((a = getPrefixAssignment()) != nullptr) {
      int size = a[0];
      assignments.push_back(Assignment {
          vector<int>(a+1,a+size+1),
            vector<int>(a+size+1,a+2*size+1),
            a[2*size+1] });
    }
  }


  
#ifdef WITH_OPEN_MPI
  void Reducer::getPrefixAssignmentsMaster(vector<Reducer::Assignment>& assignments) {
    assert(mpiRank == 0);

    // special case: empty prefix: just quit
    if (prefixSequence.size() == 0)
      return;
    
    vector<char> idleWorkers(mpiSize, 0);

    // returns the worker number and appropriate level
    auto findFirstIdleWorker = [&]()->std::pair<int,int> {
      for (size_t i = 1; i < idleWorkers.size(); ++i) {
        if (idleWorkers[i]) {
          if (lvlsForSlaves.empty()) {
            // the simple case
            // return a non-negative dummy level
            return std::make_pair(work->empty() ? -1 : i, 1);
          }
          else {
            // hierarchical case
            for (int lvl : lvlsForSlaves[i]) {
              if (!work->empty(lvl))
                return std::make_pair(i,lvl);
            }
          }
        }
      }

      return std::make_pair(-1,-1);
    };

    auto allSlavesIdle = [&idleWorkers]()->bool {
      bool b = true;
      for (size_t i = 1; i < idleWorkers.size(); ++i)
        b = b && idleWorkers[i];
      return b;
    };

    // initialize stack
    initWork();

    // check whether we should pop and update firstIdleWorker and firstIdleWorkerLvl
    std::pair<int,int> firstIdleWorkerAndLvl = std::make_pair(-1,-1);
    int& firstIdleWorker = firstIdleWorkerAndLvl.first;
    int& firstIdleWorkerLvl = firstIdleWorkerAndLvl.second;
    auto shouldPop = [&]()->bool {
      firstIdleWorkerAndLvl = findFirstIdleWorker();
      return firstIdleWorker > 0 && firstIdleWorkerLvl >= 0;
    };

      
    while (true) {
      // check if we are done
      if (work->empty() && allSlavesIdle()) {
        for (size_t i = 1; i < idleWorkers.size(); ++i)
          mpi::sendStatus(i, mpi::PLEASE_QUIT, mpiRank);
        break;
      }
      
      while (shouldPop()) {
        // we should pop something to the first idle worker
        StackElement e;
        assert(work->pop(e, firstIdleWorkerLvl));
        mpi::sendInts(e.serialize(), firstIdleWorker, mpi::PLEASE_WORK, mpiRank);
        idleWorkers[firstIdleWorker] = false;
      }

      MPI_Status status;
      mpi::recvInts(mpiBuffer, MPI_ANY_SOURCE, MPI_ANY_TAG, &status, mpiRank);
      
      if (status.MPI_TAG == mpi::PLEASE_PUSH) {
        work->push(StackElement::deserialize(mpiBuffer));
      }
      else if (status.MPI_TAG == mpi::PLEASE_MARK_ME_IDLE) {
        idleWorkers[status.MPI_SOURCE] = 1;
      }
      else if (status.MPI_TAG == mpi::PLEASE_STORE_ASSIGNMENT) {
        int* a = &mpiBuffer[0];
        int size = a[0];
        assignments.push_back(Assignment {
            vector<int>(a+1,a+size+1),
              vector<int>(a+size+1,a+2*size+1),
              a[2*size+1] });
      }
      else {
        throw runtime_error("MPI master " + to_string(mpiRank) + ": Invalid tag encountered: " +
                            mpi::tagToString(status.MPI_TAG));
      }
    }  
  }

  void Reducer::getPrefixAssignmentsSlave() {
    const int* a = nullptr;
    while ((a = getPrefixAssignment()) != nullptr) {
      mpi::sendInts(a, 2*a[0]+2, 0, mpi::PLEASE_STORE_ASSIGNMENT, mpiRank);
    }
  }
#endif // WITH_OPEN_MPI

  

  void Reducer::initWork() {
    if (prefixSequence.size() == 0)
      return; // nothing to do

    /* Initialize minimum indicators for orbits of the base graph. */
    // note: only do this if we have non-empty prefix sequence
    vector<int> seedMin = orbitMinimumIndicators(base, nullptr);
    
    int p = 0;
    int firstVar = -1;
    for(; p < static_cast<int>(traversals[0].size()); p++) {
      if (seedMin[traversals[0][p](prefixSequence[0])]) {
        firstVar = traversals[0][p](prefixSequence[0]);
        break;
      }
    }
    if (p == static_cast<int>(traversals[0].size()))
      throw runtime_error("no minimum found for base orbit");
    assert(p == 0); // this need not be the case but it seems that there are
    // ''unexpected guarantees'' with respect to the
    // construction of traversals and the construction of the
    // orbits that actually guarantee this?
    assert(firstVar != -1);

    vector<int> vars { firstVar };
    vector<int> vals { 0 };
    pushNext(&seedMin[0], &vars[0], &vals[0], 1);
  }
}

