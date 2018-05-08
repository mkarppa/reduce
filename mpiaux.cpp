#include "mpiaux.hpp"
#include <stdexcept>
#include <iostream>

using std::vector;
using std::runtime_error;
using std::to_string;
using std::string;
using std::endl;

namespace reduce {
#ifdef WITH_OPEN_MPI
  namespace mpi {
    static size_t mpiSentMessagesCount = 0,
      mpiReceivedMessagesCount = 0,
      mpiSentBytesCount = 0,
      mpiReceivedBytesCount = 0;
    static std::chrono::nanoseconds::rep totalCommunicationTime = 0;


    
    std::string intToMpiErrorString(int err) {
      switch (err) {
      case MPI_ERR_COMM:
        return "MPI_ERR_COMM";
        break;
      case MPI_ERR_COUNT:
        return "MPI_ERR_COUNT";
        break;
      case MPI_ERR_TYPE:
        return "MPI_ERR_TYPE";
        break;
      case MPI_ERR_TAG:
        return "MPI_ERR_TAG";
        break;
      case MPI_ERR_RANK:
        return "MPI_ERR_RANK";
        break;
      default:
        return "UNKNOWN ERROR";
      }
    }


    
    static void throwMpiError(int err, int mpiRank) {
      throw runtime_error("MPI node " + to_string(mpiRank) + ": Error " +
                          intToMpiErrorString(err) + " detected!");
    }

    
    
    void sendInts(const int* s, int count, int dst, int tag, int mpiRank) {
      auto beforeSend = std::chrono::steady_clock::now();
      int err = MPI_Send(s, count, MPI_INT, dst, tag, MPI_COMM_WORLD);
      auto afterSend = std::chrono::steady_clock::now();
      std::chrono::nanoseconds diff = afterSend - beforeSend;
      totalCommunicationTime += diff.count();
      if (err != MPI_SUCCESS)
        throwMpiError(err, mpiRank);
      ++mpiSentMessagesCount;
      mpiSentBytesCount += sizeof(int)*count;

    }

    

    void sendInts(const std::vector<int>& v, int dst, int tag, int mpiRank) {
      sendInts(&v[0], v.size(), dst, tag, mpiRank);
    }



    void sendStatus(int dst, int tag, int mpiRank) {
      sendInts(nullptr, 0, dst, tag, mpiRank);
    }


    
    void recvInts(std::vector<int>& s, int src, int tag,
                     MPI_Status* status, int mpiRank) {
      auto beforeRecv = std::chrono::steady_clock::now();
      int err = MPI_Recv(&s[0], s.size(), MPI_INT, src, tag,
                         MPI_COMM_WORLD, status);
      auto afterRecv = std::chrono::steady_clock::now();
      std::chrono::nanoseconds diff = afterRecv - beforeRecv;
      totalCommunicationTime += diff.count();
      if (err != MPI_SUCCESS)
        throwMpiError(err, mpiRank);
      int count;
      err = MPI_Get_count(status, MPI_INT, &count);
      if (err != MPI_SUCCESS)
        throwMpiError(err, mpiRank);    
      ++mpiReceivedMessagesCount;
      mpiReceivedBytesCount += sizeof(int)*count;
    }

    

    string tagToString(int tag) {
      switch(tag) {
      case mpi::PLEASE_PUSH:
        return "PLEASE_PUSH";
        break;
      case mpi::PLEASE_STORE_ASSIGNMENT:
        return "PLEASE_STORE_ASSIGNMENT";
        break;
      case mpi::PLEASE_QUIT:
        return "PLEASE_QUIT";
        break;
      case mpi::PLEASE_WORK:
        return "PLEASE_WORK";
        break;
      case mpi::PLEASE_MARK_ME_IDLE:
        return "PLEASE_MARK_ME_IDLE";
        break;
      default:
        return "INVALID TAG";
        break;
      }
    }


    
    void printMpiStatistics(std::ostream& out) {
      int rank, size;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      
      out << "MPI rank: " << rank << endl;
      out << "MPI size: " << size << endl;
      out << "MPI sent messages count: " << mpiSentMessagesCount << endl;
      out << "MPI received messages count: " << mpiReceivedMessagesCount
          << endl;
      out << "MPI sent bytes count: " << mpiSentBytesCount << endl;
      out << "MPI received bytes count: " << " bytes received" << endl;
      out << "MPI total communication time: " << totalCommunicationTime
          << " ns" << std::endl;
    }    
  }
#endif // WITH_OPEN_MPI
}
