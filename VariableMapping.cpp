#include "VariableMapping.hpp"
#include "common.hpp"
#include <regex>

using std::regex;
using std::regex_match;
using std::runtime_error;
using std::string;
using std::smatch;
using std::vector;
using std::to_string;

namespace reduce {
  const Graph::Vertex VariableMapping::NOT_A_VARIABLE = 0xffffffff; 
  
  void VariableMapping::parse(std::istream& in, int n) {
    //
    // Parse variables and values
    //
   
    regex e("p variable ([0-9]+)");
    string line;
    smatch sm;
    
    if (!getline(in, line) ||
        !regex_match(line,sm,e))
      throw runtime_error("parse error -- variable format line expected");
    int v = std::stoi(sm[1]);
    if(v < 1)
      throw runtime_error("bad variable parameter v = " +
                          std::to_string(v));
    var = vector<Graph::Vertex>(v);
    varLegend = vector<string>(v);

    e = regex("v ([0-9]+) (.+)");
    for (int i = 0; i < v; i++) {
      if (!getline(in,line) ||
          !regex_match(line,sm,e))
        throw std::runtime_error("parse error -- variable line expected");
      Graph::Vertex u = stoul(sm[1]);
      if(u < 1 || u > n)
        throw runtime_error("bad variable identifier u = " + std::to_string(u));
      var[i] = u-1;
      varLegend[i] = sm.str(2);
    }

    e = regex("p value ([0-9]+)");
    if (!getline(in,line) ||
        !regex_match(line,sm,e))
      throw runtime_error("parse error -- value format line expected");
    int d = stoi(sm[1]); 
    if(d < 1)
      throw runtime_error("bad value parameter r = " + to_string(d));
    val = vector<Graph::Vertex>(d);
    valLegend = vector<string>(d);
   
    e = regex("r ([0-9]+) (.+)");
    for(int i = 0; i < d; i++) {
      if (!getline(in,line) ||
          !regex_match(line,sm,e))
        throw runtime_error("parse error -- value line expected");
      Graph::Vertex u = stoul(sm[1]);
      if(u < 1 || u > n)
        throw runtime_error("bad value identifier u = " + to_string(u));
      val[i] = u-1;
      valLegend[i] = sm.str(2);
    }
  }



  void VariableMapping::fromCnf(const CNF& cnf) {
    int v = cnf.getVariableCount();
    var = vector<Graph::Vertex>(v);
    varLegend = vector<string>(v);
    for(int i = 0; i < v; ++i) {
      var[i] = i;
      varLegend[i] = to_string(i+1);
    }

    val = vector<Graph::Vertex> { 3*v, 3*v+1 };
    valLegend = vector<string> { "false", "true" }; 
  }

  
  void VariableMapping::buildTrans(const CNF& cnf, const Graph& graph) {
    varTrans = vector<CNF::Variable>(graph.getOrder());
    if(!cnf.empty()) {
      /* Build the translation array from graph variable vertices
       * to selected CNF variables. */

      vector<Graph::Vertex> q(getVariableCount(),NOT_A_VARIABLE);
      for(size_t i = 0; i < getVariableCount(); ++i) {
        size_t u = stoul(varLegend[i])-1;

        if (u >= getVariableCount())
          throw runtime_error("parsed CNF variable in legend (" + to_string(u + 1) +
                              ") is out of range");
        varTrans[var[i]] = u;
        q[i] = u;
      }

      std::sort(q.begin(), q.end());
      for(size_t i = 1; i < getVariableCount(); i++)
        if(q[i-1] == q[i])
          throw runtime_error("repeated CNF variable (" + to_string(q[i]+1) +
                              ") in legend");
                
      /* Build the translation array from graph false/true vertices
       * to CNF values, i.e. make sure false and true are present
       * and in this order. */

      if (getValueCount() != 2)
        throw runtime_error("value range does not consist of 'false' and "
                            "'true'");
      
      if (valLegend[0] == "true" && valLegend[1] == "false") {
        std::swap(val[0],val[1]);
        std::swap(valLegend[0],valLegend[1]);
      }
      else if (valLegend[1] != "true" ||
               valLegend[0] != "false") {
        throw runtime_error("value range does not consist of 'false' and 'true'");
      }
    }
    else {
      for (size_t i = 0; i < getVariableCount(); ++i)
        varTrans[i] = i;
    }


      /* Test variables for repeated elements. */
    vector<Graph::Vertex> q(var.begin(), var.end());
    std::sort(q.begin(), q.end());
    for(size_t i = 1; i < q.size(); ++i)
      if(q[i-1] == q[i])
        throw runtime_error("variable list repeats an element (" +
                            to_string(q[i]+1) + ")");
    /* Test values for repeated elements. */
    q = vector<Graph::Vertex>(val.begin(), val.end());
    std::sort(q.begin(), q.end());
    for(size_t i = 1; i < q.size(); ++i)
      if(q[i-1] == q[i])
        throw runtime_error("value list repeats an element (" +
                            to_string(q[i] + 1) + ")");

  }
}
