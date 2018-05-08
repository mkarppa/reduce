#include "Permutation.hpp"
#include <stdexcept>
#include <set>
#include <string>

using std::vector;
using std::runtime_error;
using std::to_string;
using std::set;

namespace reduce {
  static bool isPermutation(const vector<int>& p) {
    set<int> q(p.begin(), p.end());
    for (size_t i = 0; i < p.size(); ++i) 
      if (!q.count(i))
        return false;

    return true;   
  }


  
  Permutation::Permutation(size_t n) : Permutation([n]()->vector<int> {
      vector<int> q(n);
      for (size_t i = 0; i < n; ++i)
        q[i] = i;
      return q;
    }()) {
  }

  

  Permutation::Permutation(const std::vector<int>& p) : p(p) {
    if (!isPermutation(p))
      throw runtime_error("Invalid permutation");
  }

  int Permutation::operator()(int v) const {
    if (v >= static_cast<int>(p.size()) || v < 0)
      throw runtime_error("Invalid argument");
    return p[v];
  }

  std::vector<int> Permutation::operator()(const std::vector<int>& v) const {
    if (v.size() != p.size())
      throw runtime_error("invalid argument");
    vector<int> w(v.size());
    for (size_t i = 0; i < p.size(); ++i) {
      if (v[i] < 0 || v[i] >= static_cast<int>(p.size()))
        throw runtime_error("invalid argument");
      w[i] = p[v[i]];
    }
    return w;
  }

  Permutation Permutation::operator*(const Permutation& that) const {
    if (p.size() != that.p.size())
      throw runtime_error("Invalid argument");
    vector<int> q(p.size());
    for (size_t i = 0; i < p.size(); ++i)
      q[i] = p[that.p[i]];
    return Permutation(q);
  }
}

