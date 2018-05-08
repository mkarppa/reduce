#ifndef REDUCE2_VARIABLEMAPPING_HPP
#define REDUCE2_VARIABLEMAPPING_HPP

#include "Graph.hpp"


namespace reduce { 
  class VariableMapping {
  public:
    // constant used to signal that a vertex is not a variable vertex
    static const Graph::Vertex NOT_A_VARIABLE;

    // create an empty mapping
    VariableMapping() = default;
    
    // parse from input and build translation array to the selected CNF
    // variables
    explicit VariableMapping(std::istream& in, const CNF& cnf, const Graph& graph) {
      parse(in, graph.getOrder());
      buildTrans(cnf, graph);
    }

    // build directly from CNF
    explicit VariableMapping(const CNF& cnf, const Graph& graph) {
      fromCnf(cnf);
      buildTrans(cnf, graph);
    }

    size_t getVariableCount() const {
      return var.size();
    }

    size_t getValueCount() const {
      return val.size();
    }

    // returns the ith variable vertex
    Graph::Vertex getVariable(size_t i) const {
      return var[i];
    }

    // returns the ith value vertex
    Graph::Vertex getValue(size_t i) const {
      return val[i];
    }

    // translates a variable from graph to CNF
    inline CNF::Variable graphToCnf(Graph::Vertex u) const {
      return varTrans[u];
    }

    // Return the string identifier for the given variable
    inline const std::string& variableLegend(size_t v) const {
      return varLegend[v];
    }

    // Return the string identifier for the given value
    inline const std::string& valueLegend(size_t v) const {
      return valLegend[v];
    }

    
  private:
    void fromCnf(const CNF& cnf);
    // parse from dimacs-like input
    // n = order of the base graph
    void parse(std::istream& in, int n);
    void buildTrans(const CNF& cnf, const Graph& graph);

    std::vector<Graph::Vertex> var;       // Variable vertices in base graph
    std::vector<Graph::Vertex> val;       // Value vertices in base graph
    std::vector<std::string>   varLegend; // String identifiers for the variables
    std::vector<std::string>   valLegend; // String identifiers for the values
    std::vector<CNF::Variable> varTrans;  // Translation from graph to CNF variables

  };
}

#endif // REDUCE2_VARIABLEMAPPING_HPP
