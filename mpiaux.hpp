#ifndef REDUCE2_MPIAUX_HPP
#define REDUCE2_MPIAUX_HPP

#ifdef WITH_OPEN_MPI
#include <mpi.h>
#endif // WITH_OPEN_MPI

#include <vector>
#include <string>
#include <iostream>
#include <chrono>


namespace reduce {
#ifdef WITH_OPEN_MPI
  // MPI auxiliary functions
  namespace mpi {
    void sendInts(const int* s, int count, int dst, int tag, int mpiRank);
    void sendInts(const std::vector<int>& v, int dst, int tag, int mpiRank);
    void sendStatus(int dst, int tag, int mpiRank);
    void recvInts(std::vector<int>& s, int src, int tag,
                  MPI_Status* status, int mpiRank);
    std::string intToMpiErrorString(int err);
    std::string tagToString(int tag);
    void printMpiStatistics(std::ostream&);

    // message tags
    const int PLEASE_PUSH             = 1;
    const int PLEASE_QUIT             = 2;
    const int PLEASE_WORK             = 3;
    const int PLEASE_MARK_ME_IDLE     = 4;
    const int PLEASE_STORE_ASSIGNMENT = 5;
  }
#endif // WITH_OPEN_MPI
}

#endif // REDUCE2_MPIAUX_HPP
