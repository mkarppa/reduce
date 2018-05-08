#ifndef REDUCE2_PERMUTATION_HPP
#define REDUCE2_PERMUTATION_HPP

#include <vector>
#include <cstdlib>

namespace reduce {
  class Permutation {
  public:
    // construct an identity permutation on n vertices
    explicit Permutation(size_t n);

    // construct a custom permutation on vertices (inferred from the size of
    // the vector)
    explicit Permutation(const std::vector<int>& p);

    explicit Permutation() = default;

    // permute a single element
    int operator()(int v) const;

    // permute a vector
    std::vector<int> operator()(const std::vector<int>& v) const;

    // left-hand composition
    Permutation operator*(const Permutation& that) const;

    // use this only for debuggery, ordinarily the application operators ought
    // to be used
    const std::vector<int>& getP() const {
      return p;
    }
    
  private:
    std::vector<int> p;
  };
}

#endif // REDUCE2_PERMUTATION_HPP
