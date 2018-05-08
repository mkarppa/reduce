CXX=mpic++
#CXX=c++

NAUTY_PATH=./nauty/nauty26r10

NAUTY_A=$(NAUTY_PATH)/nauty.a

GMP_PATH=./gmp/gmp-6.1.2

GMP_A=$(GMP_PATH)/.libs/libgmp.a

TCLAP_PATH=./tclap/tclap-1.2.1/include

all: reduce

CXXFLAGS =-O3 -std=c++14 -Wall -pedantic -Wextra -I$(NAUTY_PATH) -I$(GMP_PATH) \
	-I$(TCLAP_PATH) -c -DCOMMITID=\"$(COMMITID)\" -DWITH_OPEN_MPI

COMMITID=$(shell git rev-parse HEAD)

HDRS=Graph.hpp common.hpp Reducer.hpp CNF.hpp VariableMapping.hpp Permutation.hpp Stack.hpp mpiaux.hpp

OBJS=reduce.o Graph.o common.o Reducer.o CNF.o VariableMapping.o \
	Permutation.o ConflictChecker.o Stack.o mpiaux.o

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $<

reduce: $(OBJS)
	$(CXX) $(LDFLAGS) -o reduce $(OBJS) \
	$(GMP_A) \
	$(NAUTY_A)

clean:
	rm -f reduce *.o *~ *.log 
