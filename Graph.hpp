/********************************************************* Graph operations. */

#ifndef REDUCE2_GRAPH_HPP
#define REDUCE2_GRAPH_HPP

#include "CNF.hpp"
#include "Permutation.hpp"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <gmpxx.h>
extern "C" {
#include "nausparse.h"
}
#include <chrono>

/******************************************************* External interface. */

namespace reduce {
  // forward declarations
  class VariableMapping;
  
  /* Graph data type. */
  class Graph {
  public:
    typedef int32_t Vertex;
    // two vertices concatenated, numerically smaller in the upper bits
    typedef int64_t Edge; 
    typedef int Size; // record sizes as int (for compatibility with nauty)


    /**
     * Parse graph from input stream
     */
    explicit Graph(std::istream& in) {
      parse(in);
    }



    /**
     * Create an empty graph of given order
     */
    explicit Graph(Size order) {
      init(order);
    }



    /**
     * Construct graph from CNF
     */
    explicit Graph(const CNF& cnf);



    /**
     * Construct an empty graph
     */
    explicit Graph() :
      Graph(0) {
    }


    
    ~Graph() = default;
    Graph(const Graph&) = default;
    Graph& operator=(const Graph&) = default;
    Graph(Graph&& that) = default;
    Graph& operator=(Graph&& that) = default;


    
    /********************************************* Returns the order of a graph. */
    Size getOrder() const {
      return order;
    }

    

    Size getNumEdges() const {
      return edgeBuf.size();
    }


    
    // relabels the graph
    Graph relabel(const std::vector<int>& p);
    Graph relabelInv(const std::vector<int>& p);


    
    // add an edge to the graph
    void addEdge(Vertex i, Vertex j);


    
    /*************************************** Computes a graph in canonical form. */
    Graph canForm();

    

    void printOrbits(std::ostream& out, const VariableMapping& varMap) const;

    

    // Returns automorphism orbits
    const std::vector<int>& getOrbits() const;

    

    // returns automorphism group generators
    const std::vector<std::vector<int> >& getAutGen() const;

    const std::vector<Permutation>& getAutGenPerm() const;

    

    // computes the size of the automorphism group (a big integer)
    mpz_class getAutOrder() const;

    

    /************ Tests whether two vertices are in the same automorphism orbit. */
    bool sameOrbit(Vertex i, Vertex j) const;

    

    // Build the orbit cells
    const std::vector<int>& orbitCells() const;
    

    
    /***************** Returns a stabilizer sequence for the automorphism group. */
    const std::vector<int>& getStabSeq() const;

   

    // Computes canonical labeling for a graph
    const std::vector<int>& canLab();    



    // computes the orbit of a given vertex
    std::vector<Vertex> orbitOf(Vertex v) const;
    
  private:
    Size order; // number of vertices
    std::vector<Edge> edgeBuf;    // edge buffer
    mutable std::vector<Edge> canEdgeBuf; // canonical edges
    bool edgeBufIsSorted;    // whether the edge buffer has been sorted
    mutable std::vector<int> orb; // automorpihsm orbits
    mutable std::vector<int> lab; // vertex labels (vertices in "some order", see nauty docs)
    mutable std::vector<int> ptn; // vertex partitions (binary array where 1 indicates the
                                  //                    end of color, see nauty docs)
    mutable std::vector<int> orbCells; // orbit cells
    mutable std::vector<std::vector<int> > autGen; // automorphism group generators
    mutable std::vector<Permutation> autGenPerm;
    mutable std::vector<int> autIdx;
    mutable std::vector<int> stabSeq;
    mutable bool haveCan;          // whether the canonical form has been computed

    // A simple graph parser
    void parse(std::istream& in);

    // initialize the data structures
    void init(Size order);

    // compute the canonical form if not already computed
    void getCan() const;
    
    // swap the two graphs
    void internalSwap(Graph& that);

    void sortEdgeBuf();

    // nauty user-defined procedures below
    static void lvlproc(int *lab, int *ptn, int lvl, int *orb, statsblk *stats,
                        int tv, int idx, int tcellsize, int numcells,
                        int childcount, int n);
    static void automproc(int numgen, int *p, int *orb, 
                          int numorb, int stabv, int n);
  };

  // returns total time spent in nauty calls in nanoseconds
  std::chrono::nanoseconds::rep getTotalNautyTime();
  size_t getNautyCallCount();
}

  
#endif // REDUCE2_GRAPH_HPP
