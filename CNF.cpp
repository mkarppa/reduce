#include "CNF.hpp"
#include <regex>
#include <sstream>

using std::string;
using std::regex;
using std::regex_match;
using std::smatch;
using std::runtime_error;
using std::to_string;
using std::vector;
using std::cerr;
using std::endl;

namespace reduce {
  // eat comments
  static std::istream& getNextNonCommentLine(std::istream& is,
                                             std::string& out) {
    while (getline(is, out)) {
      if (out.length() == 0 || out[0] != 'c')
        return is;
    }
    return is;
  }
  
  void CNF::parseDimacs(std::istream& is) {
    string line;
    size_t nv = 0, nc = 0;
    regex rheader("p cnf ([1-9][0-9]*) ([1-9][0-9]*)");
    regex rclause("(-?[1-9][0-9]* )*0");

    // read header
    if (!getNextNonCommentLine(is,line))
      throw runtime_error("Unexpected end of input!");
    
    smatch sm;
      
    if (regex_match(line,sm,rheader)) {
      nv = stoul(sm[1]);
      nc = stoul(sm[2]);
      for (size_t v = 1; v <= nv; ++v) {
        variables.insert(v);
        literalOccurrences[v]  = std::vector<ClauseIdx>();
        literalOccurrences[-v] = std::vector<ClauseIdx>();
      }
    }
    else {
      throw std::runtime_error("Invalid DIMACS encountered: "
                               "Header expected, none given.");
    }

    size_t i = 0;
    while (i < nc) {
      if (!getline(is, line)) 
        throw runtime_error("Unexpected end of input!");

      if (line.length() > 0 && line[0] == 'c') {
        continue; // ignore comments
      }
      else if (regex_match(line,rclause)) {
        std::istringstream iss(line);
        Clause cl;
        Literal v;
        while (iss >> v) {
          if (v == 0) {
            break;
          }
          else if (v > static_cast<int>(nv)) {
            throw std::runtime_error("Invalid DIMACS encountered: literal `" +
                                     to_string(v) +
                                     "' exceeds variable count of `" +
                                     to_string(nv) + "'");
          }
          else {
            cl.push_back(v);
            literalOccurrences[v].push_back(clauses.size());
          }
        }
        clauses.push_back(cl);
        ++i;
      }
      else {
        throw std::runtime_error("Invalid DIMACS encountered: "
                                 "malformed line `" + line + "'");
      }
    }
    
    if (clauses.size() != nc) {
      throw std::runtime_error("Invalid DIMACS encountered: "
                               "clause count mismatch, expected " +
                               to_string(nc) + " clauses but got " +
                               to_string(clauses.size()));
    }
  }
}
