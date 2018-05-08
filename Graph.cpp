/********************************************************* Graph operations. */

#include "Graph.hpp"
#include "VariableMapping.hpp"
#include "common.hpp"
#include <algorithm>
#include <regex>
#include <cassert>

using std::string;
using std::regex;
using std::regex_match;
using std::smatch;
using std::to_string;
using std::stoi;
using std::runtime_error;
using std::vector;
using std::swap;
using std::cerr;
using std::endl;

// we will assume this for nauty
static_assert(std::is_same<int32_t,int>::value, "expecting int32_t and int to "
              "be the same type for nauty");

namespace reduce {

  static std::chrono::nanoseconds::rep totalNautyTime = 0;
  size_t nautyCallCount = 0;
  
  std::chrono::nanoseconds::rep getTotalNautyTime() {
    return totalNautyTime;
  }

  size_t getNautyCallCount() {
    return nautyCallCount;
  }
  
  /********************************************************** Graph data type. */


  /************************************* Initialization and release functions. */
 
  void Graph::init(Size n) {
    order = n;
    edgeBuf = vector<Edge>();
    canEdgeBuf = vector<Edge>();
    edgeBufIsSorted = false;
    lab        = vector<int>(order);
    ptn        = vector<int>(order,1);
    orb        = vector<int>(order);
    for(Size i = 0; i < order; ++i) {
      lab[i] = i;
    }
    if (order > 0)
      ptn[order-1] = 0;    
    haveCan = false;
  }

 


  /*************************************** Subroutines for working with edges. */

  static Graph::Edge edge_make(Graph::Vertex i, Graph::Vertex j) {
    if(i < j)
      return ((static_cast<Graph::Edge>(i))<<32) | j;
    else 
      return ((static_cast<Graph::Edge>(j))<<32) | i;
  }

  

  static Graph::Vertex edge_i(Graph::Edge e) {
    return e >> 32; 
  }

  

  static Graph::Vertex edge_j(Graph::Edge e) {
    return e & 0xFFFFFFFF;  
  }

  

  static Graph::Edge edge_relabel(const vector<int>& p, Graph::Edge e) {
    return edge_make(p[edge_i(e)], p[edge_j(e)]);
  }


  

  void Graph::sortEdgeBuf() {
    if(!edgeBufIsSorted) {
      std::sort(edgeBuf.begin(), edgeBuf.end());
      for (Size l = 1; l < getNumEdges(); ++l) 
        if (edgeBuf[l-1] >= edgeBuf[l]) 
          throw runtime_error("found repeated edge or bad edge sort");      
    }
    edgeBufIsSorted = true;
  }



  /********************************************************** Relabel a graph. */

  // verify that we're dealing with a permutation
  static void graph_permcheck(int n, const vector<int>& p) {
    vector<int> pc(n,0);
    for(int i = 0; i < n; i++) {
      if(p[i] < 0 || p[i] >= n || pc[p[i]] != 0)
        throw runtime_error("invalid permutation");
      pc[p[i]] = 1;
    }
  }


  
  Graph Graph::relabel(const std::vector<int>& p) {
    graph_permcheck(order, p);
    Graph r(order);
    for(Size i = 0; i < order; i++) {
      r.lab[i] = p[lab[i]];
      r.ptn[i] = ptn[i];
    }
    for (Edge e : edgeBuf) {
      r.edgeBuf.push_back(edge_relabel(p,e));
    }
    r.edgeBufIsSorted = false;
    return r;
  }


  
  Graph Graph::relabelInv(const std::vector<int>& p) {
    graph_permcheck(getOrder(), p);
    vector<int> pinv(getOrder());
    for(Size i = 0; i < getOrder(); i++)
      pinv[p[i]] = i;
    return relabel(pinv);
  }


  
  /************************************************************** Add an edge. */

  void Graph::addEdge(Vertex i, Vertex j) {
    if (i >= getOrder() || j >= getOrder() || i == j) {
      throw runtime_error("bad edge (i = " + to_string(i) +
                          ", j = " + to_string(j));
    }
    haveCan = false;
    edgeBufIsSorted = false;
    edgeBuf.push_back(edge_make(i,j));
  }



  static const Graph* autom_g;

  void Graph::lvlproc(int* /* lab */,
                      int* /* ptn */,
                      int /* lvl */,
                      int* /* orb */,
                      statsblk* /* stats */,
                      int tv, int idx, int /* tcellsize */, int /* numcells */,
                      int /* childcount */, int /* n */) {
    autom_g->stabSeq.push_back(tv);
    autom_g->autIdx.push_back(idx);
  }

  

  void Graph::automproc(int /* numgen */, int* p, int* /* orb */, 
                        int /* numorb */, int /* stabv */, int /* n */) {
    autom_g->autGen.push_back(vector<int>(p, p + autom_g->getOrder()));
  }

  

  void Graph::getCan() const {
    if (haveCan)
      return;

    push_time();

    Size n = getOrder();
    Size m = getNumEdges();
    
    autGen.clear();
    autGenPerm.clear();
    autIdx.clear();

    sparsegraph ng, ncg;
    SG_INIT(ng);
    SG_INIT(ncg);
    static DEFAULTOPTIONS_SPARSEGRAPH(options);
    statsblk stats;

    int mm = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE, mm, n, NAUTYVERSIONID);

    vector<size_t> v(n);
    vector<int> d(n,0);
    vector<int> e(m*2);
    vector<size_t> cv(n); 
    vector<int> cd(n);
    vector<int> ce(m*2);

    ng.nv   = n;
    ng.nde  = m*2;
    ng.vlen = n;
    ng.dlen = n;
    ng.elen = m*2;
    ng.wlen = 0;
    ng.v    = &v[0];
    ng.d    = &d[0];
    ng.e    = &e[0];
    ng.w    = nullptr;

    ncg.nv   = n;
    ncg.nde  = m*2;
    ncg.vlen = n;
    ncg.dlen = n;
    ncg.elen = m*2;
    ncg.wlen = 0;
    ncg.v    = &cv[0];
    ncg.d    = &cd[0];
    ncg.e    = &ce[0];
    ncg.w    = nullptr;

    for (Edge b : edgeBuf) {
      Vertex i = edge_i(b);
      Vertex j = edge_j(b);
      d[i]++;
      d[j]++;
    }
    v[0] = d[0];
    for (Size i = 1; i < n; ++i)
      v[i] = v[i-1] + d[i];
    if (v[n-1] != 2*static_cast<size_t>(m))
      throw runtime_error("bad v array");
    
    for (Edge b : edgeBuf) {
      Vertex i = edge_i(b);
      Vertex j = edge_j(b);
      e[--v[j]] = i;
      e[--v[i]] = j;
    }
    
    options.defaultptn    = 0;
    options.getcanon      = 1;
    options.userautomproc = &automproc;
    options.userlevelproc = &lvlproc;

    autom_g = this;
    assert(autGen.size() == 0);

    auto beforeNauty = std::chrono::steady_clock::now();
    sparsenauty(&ng, &lab[0], &ptn[0], &orb[0], &options, &stats, &ncg);
    auto afterNauty = std::chrono::steady_clock::now();
    std::chrono::nanoseconds diff = afterNauty - beforeNauty;
    totalNautyTime += diff.count();
    ++nautyCallCount;
    
    Size l = 0;
    canEdgeBuf.clear();
    
    for(Vertex i = 0; i < n; ++i) {
      std::sort(ce.begin() + cv[i],
                ce.begin() + cv[i] + cd[i]);
      for (auto it = ce.begin() + cv[i];
           it != ce.begin() + cv[i] + cd[i];
           ++it) {
        Vertex j = *it;
        if(i < j) {
          canEdgeBuf.push_back(edge_make(i, j));
          ++l;
        }
      }
    }
    assert(canEdgeBuf.size() == static_cast<size_t>(l));
    assert(canEdgeBuf.size() == static_cast<size_t>(m));
    if (l != m)
      throw runtime_error("bad canonical form (l = " + to_string(l) +
                          ", m = " + to_string(m) + ")");
    for (l = 1; l < m; l++)
      if (canEdgeBuf[l-1] >= canEdgeBuf[l])
        throw runtime_error("bad sort in canonical form");

    autGenPerm = vector<Permutation>(autGen.size());
    for (size_t i = 0; i < autGen.size(); ++i)
      autGenPerm[i] = Permutation(autGen[i]);

    pop_print_time("nauty");

    haveCan = true;
  }

  

  const std::vector<int>& Graph::canLab() {
    getCan();
    return lab;
  }

  

  Graph Graph::canForm() {
    getCan();
    Graph cg(getOrder());
    cg.edgeBuf = canEdgeBuf;
    cg.edgeBufIsSorted = true;
    for(Size i = 0; i < getOrder(); i++)
      cg.ptn[i] = ptn[i];
    return cg;
  }

  

  const std::vector<std::vector<int> >& Graph::getAutGen() const {
    getCan();
    return autGen;
  }

  const std::vector<Permutation>& Graph::getAutGenPerm() const {
    getCan();
    return autGenPerm;
  }

  
  
  const std::vector<int>& Graph::getStabSeq() const {
    getCan();
    return stabSeq;
  }

  

  const std::vector<int>& Graph::getOrbits() const {
    getCan();
    return orb;
  }

  

  bool Graph::sameOrbit(Vertex i, Vertex j) const {
    if (i >= getOrder() || j >= getOrder())
      throw runtime_error("vertex number out of bounds (i = " + to_string(i) +
                          ", j = " + to_string(j));
    if (i == j)
      return true;
    
    getCan();

    if(orb[i] == orb[j])
      return true;
    else
      return false;
  }


  
  const vector<int>& Graph::orbitCells() const {
    getCan();
    size_t n = getOrder();
    if (orbCells.empty())
      orbCells = vector<int>(n);
    assert(orbCells.size() == n);
    
    for (size_t i = 0; i < n; ++i)
      orbCells[i] = i;
    
    assert(orb.size() == n);
    std::sort(orbCells.begin(), orbCells.end(), [&](int x, int y) {
        return orb[x] < orb[y];
      });
    
    size_t last = 0;

    while(last < n) {
      size_t i = last + 1;

      for(; i < n; ++i)
        if(orb[orbCells[i-1]] != orb[orbCells[i]])
          break;

      std::sort(orbCells.begin() + last, orbCells.begin() + i);
      last = i;
    }
    
    return orbCells;
  }


  
  void Graph::parse(std::istream& in) {
    string line;
    smatch sm;
    regex e("p edge ([0-9]+) ([0-9]+)");

    while (getline(in,line)) {
      // eat possible comments at the beginning
      if (line.length() > 0 && line[0] == 'c')
        continue;
      else {
        regex_match(line,sm,e);
        break;
      }
    }
    if (sm.empty())
      throw std::runtime_error("parse error -- graph format line expected");

    Size n = stoi(sm[1]);
    Size m = stoi(sm[2]);
    
    if(n <= 1)
      throw std::runtime_error("bad graph parameters n = " +
                               std::to_string(n) + ", m = " +
                               std::to_string(m));
    init(n);
    vector<int> colors(n,-1);

    e = regex("e ([0-9]+) ([0-9]+)");
    for(Size i = 0; i < m; ++i) {
      if (!getline(in,line) ||
          !regex_match(line,sm,e)) {
        throw std::runtime_error("parse error -- edge line expected");
      }
      Vertex u = std::stoul(sm[1]), v = std::stoul(sm[2]);
      if(u < 1 || v < 1 || u == v || u > n || v > n)
        throw std::runtime_error("bad edge in input u = " + std::to_string(u) +
                                 ", v = " + std::to_string(v));
      addEdge(u-1, v-1);
    }

    e = regex("c ([0-9]+) ([0-9]+)");
    for (Size i = 0; i < n; ++i) {
      if(!getline(in,line) ||
         !regex_match(line,sm,e))
        throw std::runtime_error("parse error -- color line expected");
      Vertex u = std::stoul(sm[1]);
      int c = std::stoi(sm[2]);
      if(u < 1 || c < 0 || u > n)
        throw std::runtime_error("bad color u = " + std::to_string(u) +
                                 ", c = " + std::to_string(c));
      colors[u-1] = c;
    }
    for (Size i = 0; i < n; i++)
      if(colors[i] == -1)
        throw std::runtime_error("vertex u = " + std::to_string(i+1) +
                                 " did not receive a color");

    // indirect sort of labels wrt. colors
    std::sort(lab.begin(),lab.begin()+n,[&](int x, int y) {
        return colors[x] < colors[y];
      });

    for (Size i = 0; i < n; ++i)
      if (i == n-1 || colors[lab[i]] != colors[lab[i+1]])
        ptn[i] = 0;
      else
        ptn[i] = 1;
  }

  

  /*********************************************************** Graph printing. */



  void Graph::printOrbits(std::ostream& out,
                          const VariableMapping& varMap) const {
    getCan();
    Size n = getOrder();
    Size l = varMap.getVariableCount();
    const vector<int>& p = orbitCells();
    const vector<int>& c = getOrbits();
    vector<int> q(n, 0);

    for (Size i = 0; i < l; i++) {
      if (varMap.getVariable(i) >= n)
        throw runtime_error("bad m array");
      q[varMap.getVariable(i)] = 1;
    }
    
    bool have_previous = false;
    for (Size s = 0; s < n; ++s) {
      Size u = s+1;
      while (u < n && c[p[s]] == c[p[u]])
        ++u;
      
      for (Size j = s+1; j < u; ++j)
        if (q[p[s]] != q[p[j]])
          throw runtime_error("bad m array -- not a union of orbits");
      
      if (q[p[s]]) {
        out << (have_previous ? " | " : "");
        printArray(out, vector<int>(p.begin()+s, p.begin()+u));
        have_previous = true;
      }
      s = u-1;
    }
  }



  Graph::Graph(const CNF& cnf) {
    /* Build the graph of symmetries from CNF. */
    Size nv = cnf.getVariableCount();
    Size nc = cnf.getClauseCount();
    Size n = 3*nv + 2 + nc;
    init(n);

    for (Size i = 0; i < nv; ++i) {
      addEdge(i, nv + i);
      addEdge(i, 2*nv + i);
    }

    for (CNF::ClauseIdx c = 0; c < nc; ++c) {
      for (CNF::Literal l : cnf[c]) {
        if (l < 0) {
          l = (-l)-1;
          if(l >= static_cast<int>(nv))
            throw std::runtime_error("negative literal (" +
                                     std::to_string(-(l + 1)) +
                                     ") out of range");
          addEdge(nv + l, 3*nv + 2 + c);
        }
        else {
          l = l-1;
          if(l >= static_cast<int>(nv))
            throw std::runtime_error("positive literal (" +
                                     std::to_string(l + 1) +
                                     ") out of range");
          addEdge(2*nv + l, 3*nv + 2 + c);
        }
      }
    }
      
    vector<int> colors(n,-1);
    for (Size i = 0; i < n; i++)
      colors[i] = -1;
    for (Size i = 0; i < nv; i++)
      colors[i] = 0;
    for (Size i = 0; i < nv; i++)
      colors[nv + i] = 1;
    for (Size i = 0; i < nv; i++)
      colors[2*nv + i] = 2;
    colors[3*nv + 0] = 3;
    colors[3*nv + 1] = 4;
    for (Size i = 0; i < nc; i++)
      colors[3*nv + 2 + i] = 5;        
    for (Size u = 0; u < n; u++)
      if(colors[u] == -1)
        throw runtime_error("vertex u = " + to_string(u) +
                            " did not receive a color");
    assert(lab.size() == static_cast<size_t>(n));
    std::sort(lab.begin(),lab.end(),[&colors](int x, int y) {
        return colors[x] < colors[y];
      });
    for (Size i = 0; i < n; ++i)
      if(i == n-1 || colors[lab[i]] != colors[lab[i+1]])
        ptn[i] = 0;
      else
        ptn[i] = 1;    
  }



  mpz_class Graph::getAutOrder() const {
    getCan();
    mpz_class autOrder(1L);
    for (auto ai : autIdx)
      autOrder *= ai;
    return autOrder;
  }



  std::vector<Graph::Vertex> Graph::orbitOf(Vertex v) const {
    if (v >= order)
      throw runtime_error("Invalid vertex given");

    getCan();
    
    vector<Vertex> orbit;
    for (Vertex u = 0; u < order; ++u)
      if (sameOrbit(u,v))
        orbit.push_back(u);

    return orbit;
  }
}
