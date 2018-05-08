#include "Stack.hpp"
#include "mpiaux.hpp"
#include <stdexcept>
#include <cassert>

using std::vector;
using std::endl;
using std::runtime_error;
using std::to_string;


namespace reduce {
  StackElement::StackElement(const std::vector<int>& inVars,
                             const std::vector<int>& inVals,
                             const std::vector<int>& inSeedMin) :
    vars(inVars),
    vals(inVals),
    seedMin(inSeedMin)
  {
    assert(inVars.size() == inVals.size());
  }



  Stack::~Stack() { }



  vector<int> StackElement::serialize() const {
    vector<int> out(2 + 2*vars.size() + seedMin.size());
    out[0] = vars.size();
    out[1 + 2*vars.size()] = seedMin.size();
    std::copy(vars.begin(), vars.end(), out.begin() + 1);
    std::copy(vals.begin(), vals.end(), out.begin() + 1 + vars.size());
    std::copy(seedMin.begin(), seedMin.end(),
              out.begin() + 2 + 2*vars.size());
    return out;
  }


  
  void SoloStack::push(const StackElement& e) {
    deque.push_back(e);
  } 


  
  bool SoloStack::pop(StackElement& out,int) {
    if (deque.empty())
      return false;
    out = deque.back();
    deque.pop_back();
    return true;
  }


  
  StackElement StackElement::deserialize(const vector<int>& v) {
    StackElement e;
    e.vars = vector<int>(v.begin() + 1, v.begin() + 1 + v[0]);
    e.vals = vector<int>(v.begin() + 1 + v[0], v.begin() + 1 + 2*v[0]);
    e.seedMin = vector<int>(v.begin() + 2 + 2*v[0],
                            v.begin() + 2 + 2*v[0] + v[2*v[0]+1]);
    return e;
  }

  

  bool SoloStack::empty(int) const {
    return deque.empty();
  }


  
#ifdef WITH_OPEN_MPI
  VerySimpleSlaveStack::VerySimpleSlaveStack(size_t bufferSize) :
    mpiRank(-1),
    mpiBuffer(bufferSize) {
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  }


  
  void VerySimpleSlaveStack::push(const StackElement& e) {
    const auto& s = e.serialize();
    mpi::sendInts(s,0,mpi::PLEASE_PUSH,mpiRank);
  } 


  
  bool VerySimpleSlaveStack::pop(StackElement& e,int) {
    // announce idle status, then receive answer     
    mpi::sendStatus(0, mpi::PLEASE_MARK_ME_IDLE,mpiRank);
    MPI_Status status;
    mpi::recvInts(mpiBuffer, 0, MPI_ANY_TAG, &status, mpiRank);

    if (status.MPI_TAG == mpi::PLEASE_WORK) {
      e = StackElement::deserialize(mpiBuffer);
      return true;
    }
    else if (status.MPI_TAG == mpi::PLEASE_QUIT) {
      return false;
    }
    else {
      throw runtime_error("MPI slave " + to_string(mpiRank) +
                          ": invalid MPI tag encountered: " +
                          mpi::tagToString(status.MPI_TAG));
      return false;
    }
  }


  
  bool VerySimpleSlaveStack::empty(int) const {
    return true;
  }


  
  HierarchicalSlaveStack::HierarchicalSlaveStack(const std::vector<int>& inLvls,
                                                 size_t nLvls,
                                                 size_t bufferSize) :
    lvls(nLvls,0), verySimpleMpiStack(bufferSize) {
    for (int i : inLvls)
      lvls[i] = 1;
  }


  
  HierarchicalMasterStack::
  HierarchicalMasterStack(const vector<vector<int> >& slavesForLvl) :
    deques(slavesForLvl.size()),
    slavesForLvl(slavesForLvl) {
  }


  
  bool HierarchicalMasterStack::empty(int idx) const {
    if (idx < 0) {
      for (auto& d : deques)
        if (!d.empty())
          return false;
      return true;
    }
    else
      return deques[idx].empty();
  }


  
  void HierarchicalMasterStack::push(const StackElement& e) {
    // determine level and store in appropriate deque
    int lvl = e.getSize() - 1;
    deques[lvl].push_back(e);
  }


  
  bool HierarchicalMasterStack::pop(StackElement& out, int idx) {
    // pop element from the associated deque
    if (deques[idx].empty())
      return false;

    out = deques[idx].back();
    deques[idx].pop_back();
    return true;
  }


  
  bool HierarchicalSlaveStack::empty(int idx) const {
    return solo.empty(idx);
  }


  
  void HierarchicalSlaveStack::push(const StackElement& e) {
    // check if the value should be stored locally or remotely
    int lvl = e.getSize()-1;
    if (lvls[lvl])
      solo.push(e);
    else
      verySimpleMpiStack.push(e);
  }


  
  bool HierarchicalSlaveStack::pop(StackElement& out, int idx) {
    if (!solo.empty(idx)) {
      return solo.pop(out, idx);
    }
    else {
      return verySimpleMpiStack.pop(out, idx);
    }
  }
#endif // WITH_OPEN_MPI
}

