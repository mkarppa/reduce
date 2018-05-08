#ifndef REDUCE2_STACK_HPP
#define REDUCE2_STACK_HPP

#include <vector>
#include <cstdlib>
#include <deque>
#include <iostream>

namespace reduce {
  class StackElement {
  public:
    inline int* getVars() {
      return &vars[0];
    }
    
    inline int* getVals() {
      return &vals[0];
    }

    inline int* getSeedMin() {
      return &seedMin[0];
    }
    
    inline int getSize() const {
      return vars.size();
    }

    explicit StackElement() = default;
    
    explicit StackElement(const std::vector<int>& inVars,
                          const std::vector<int>& inVals,
                          const std::vector<int>& seedMin);

    std::vector<int> serialize() const;

    static StackElement deserialize(const std::vector<int>&);

    void print(std::ostream&) const;

  private:
    std::vector<int> vars;
    std::vector<int> vals;
    std::vector<int> seedMin;

    // data layout:
    // [n][n ints][n ints][m][m ints]
    //       vars    vals     seedmin
  };

  class Stack {
  public:
    // return if the stack is empty
    // optionally request for a given substack
    virtual bool empty(int idx = -1) const = 0;
    
    virtual void push(const StackElement& e) = 0;

    // pop (optionally from a substack idx)
    virtual bool pop(StackElement&, int idx = -1) = 0;

    virtual ~Stack() = 0;
  };

  

  // stack without any MPI communication
  class SoloStack : public Stack {
  public:
    virtual bool empty(int) const override;
    
    virtual void push(const StackElement& e) override;

    virtual bool pop(StackElement&,int) override;

    virtual ~SoloStack() { }
    
  private:
    std::deque<StackElement> deque;
  };



  // stack with trivial MPI communication
  class VerySimpleSlaveStack : public Stack {
  public:
    explicit VerySimpleSlaveStack(size_t bufferSize);
    
    virtual bool empty(int) const override;
    
    virtual void push(const StackElement& e) override;

    virtual bool pop(StackElement&,int) override;

    virtual ~VerySimpleSlaveStack() { }
    
  private:
    int mpiRank;
    std::vector<int> mpiBuffer;
  }; 


  
  // stack with hierarchical MPI communication
  class HierarchicalSlaveStack : public Stack {
  public:
    explicit HierarchicalSlaveStack(const std::vector<int>& lvls,
                                    size_t nLvls,
                                    size_t bufferSize);
    
    virtual bool empty(int) const override;
    
    virtual void push(const StackElement& e) override;

    virtual bool pop(StackElement&,int) override;

    virtual ~HierarchicalSlaveStack() { }
    
  private:
    std::vector<int> lvls; // indicators of which lvls belong to this stack
    SoloStack solo;
    VerySimpleSlaveStack verySimpleMpiStack;
  }; 

  class HierarchicalMasterStack : public Stack {
  public:
    explicit HierarchicalMasterStack(const std::vector<std::vector<int> >&
                                     slavesForLvl);
    
    virtual bool empty(int) const override;
    
    virtual void push(const StackElement& e) override;

    virtual bool pop(StackElement&,int) override;

    virtual ~HierarchicalMasterStack() { }
    
  private:
    std::vector<std::deque<StackElement> > deques;
    std::vector<std::vector<int> > slavesForLvl;
  }; 
}

#endif // REDUCE2_STACK_HPP
