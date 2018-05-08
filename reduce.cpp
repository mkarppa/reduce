
/******************************************************* Symmetry reduction. */

#include "mpiaux.hpp"
#include "CNF.hpp"
#include "Reducer.hpp"
#include "common.hpp"
#include <tclap/CmdLine.h>
#include <regex>
#include <fstream>
#include <cstring>
#include <cassert>

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::to_string;

namespace {
  struct Args {
    reduce::CNF inCnf;
    reduce::Graph graph;
    reduce::VariableMapping varMap;
    vector<int> prefix;
    long threshold;
    bool symmetryOnly;
    bool incremental;
    bool verbose;
    bool printAutSizes;
    bool checkConflicts;
    string logFilename;
    vector<vector<int> > slavesForLvl;
    vector<vector<int> > lvlsForSlave;
    reduce::Stack* workStack;
    vector<double> estimationProbabilities;
  };
}

#ifdef WITH_OPEN_MPI
static vector<string> tokenize(const string& s,
                               const string& regex) {
  std::regex e(regex);
  vector<string> tokens;
  std::sregex_token_iterator begin(s.begin(), s.end(), e), end;
  std::copy(begin, end, std::back_inserter(tokens));
  return tokens;
}
#endif // WITH_OPEN_MPI

static Args parseArgs(int argc, char** argv) {   
  TCLAP::CmdLine cmd("");
  
  TCLAP::SwitchArg cnfSwitch("n","no-cnf","do not expect CNF in input",true);
  TCLAP::SwitchArg graphSwitch("g","graph","separate symmetry graph supplied in input");
  TCLAP::ValueArg<string> prefixValue("p","prefix","use the prefix <SEQ> of variable vertices",false,"","SEQ");
  TCLAP::ValueArg<long> thresholdValue("t","threshold","output partial assignment when |Aut| <= <N>",false,-1,"N");
#ifdef WITH_OPEN_MPI
  TCLAP::ValueArg<string> mpiJobPartitionValue("m","mpi-partition","MPI work partition list",false,"","string");
#endif // WITH_OPEN_MPI
  TCLAP::ValueArg<string> logFileValue("l","log-file","Statistics log filename"
#ifdef WITH_OPEN_MPI
                                       ", use `%r' as placeholder MPI rank"
#endif // WITH_OPEN_MPI
                                       ", or leave unspecified for stderr", 
                                       false,"","string");
  TCLAP::ValueArg<string> estimationValue("e","estimation-probabilities",
                                          "A string of k comma-separated "
                                          "floating point values in [0,1] for "
                                          "determining probability of "
                                          "including a branch in tree-size "
                                          "estimate", false, "", "string");
  TCLAP::SwitchArg symmetryOnlySwitch("s","symmetry-only","print symmetry information only");
  TCLAP::SwitchArg incrementalSwitch("i","incremental","give output in icnf format");
  TCLAP::SwitchArg verboseSwitch("v","verbose","verbose output");
  TCLAP::SwitchArg printAutSizesSwitch("a","print-aut-sizes","print automorphism group sizes wrt. cubes");
  TCLAP::SwitchArg conflictSwitch("c","enable-conflict-checking","enable conflict checking");
  TCLAP::UnlabeledValueArg<string> inputFilenameValue("file","input "
						      "filename or `-' for stdin",
						      true, "-", "file");
  cmd.add(cnfSwitch);
  cmd.add(graphSwitch);
  cmd.add(prefixValue);
  cmd.add(thresholdValue);
#ifdef WITH_OPEN_MPI
  cmd.add(mpiJobPartitionValue);
#endif // WITH_OPEN_MPI
  cmd.add(logFileValue);
  cmd.add(estimationValue);
  cmd.add(symmetryOnlySwitch);
  cmd.add(incrementalSwitch);
  cmd.add(verboseSwitch);
  cmd.add(printAutSizesSwitch);
  cmd.add(conflictSwitch);
  cmd.add(inputFilenameValue);
  
  cmd.parse(argc,argv);

  string infileName = inputFilenameValue.getValue();
  std::ifstream ifs;
  if (infileName != "-") {
    ifs.open(infileName);
    if (!ifs)
      throw std::runtime_error("Failed to open `" + infileName +
			       "' for read access");
  }
  std::istream& in = (infileName == "-") ? std::cin : ifs;

  Args args;

  // parse CNF if expected
  if (cnfSwitch.getValue()) {
    args.inCnf = reduce::CNF(in);
  }

  // likewise, parse graph if expected
  // also parse variable mappings
  if (graphSwitch.getValue()) {
    args.graph = reduce::Graph(in);
    args.varMap = reduce::VariableMapping(in, args.inCnf,
                                          args.graph);
  }
  else if (cnfSwitch.getValue()) {
    args.graph = reduce::Graph(args.inCnf);
    args.varMap = reduce::VariableMapping(args.inCnf, args.graph);
  }
  else
    throw TCLAP::ArgException("cannot build the symmetry graph since no CNF was given");


  if (prefixValue.getValue() != "") {
    string s = prefixValue.getValue();
    std::regex e("([1-9][0-9]* )*[1-9][0-9]*");
    if (!std::regex_match(s,e))
      throw TCLAP::ArgException("given prefix \"" + s + "\" is not a list of "
                                "positive space-separated integers", "ERROR");
    std::istringstream iss(s);
    int tmp;
    while (iss >> tmp)
      args.prefix.push_back(tmp-1);
  }

#ifdef WITH_OPEN_MPI
  int mpiSize;
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  int mpiRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  size_t prefixSize = args.prefix.size();
  size_t graphOrder = args.graph.getOrder();
  size_t mpiBufferSize = (2*prefixSize+2+graphOrder)*(prefixSize+1);
  
  if (mpiJobPartitionValue.getValue() != "") {
    args.slavesForLvl = vector<vector<int> >(args.prefix.size());
    args.lvlsForSlave = vector<vector<int> >(mpiSize);
    auto& slavesForLvl = args.slavesForLvl;
    auto& lvlsForSlave = args.lvlsForSlave;
    
    string s = mpiJobPartitionValue.getValue();
    vector<string> tokens = tokenize(s,"[0-9]+(-[0-9]+)?:[0-9]+(-[0-9]+)?");

    std::regex e("([0-9]+)-?([0-9]+)?:([0-9]+)-?([0-9]+)?");
    for (auto s : tokens) {
      std::smatch m;
      std::regex_match(s, m, e);
      if (!std::regex_match(s, m, e))
        throw TCLAP::ArgException("given prefix MPI partition string is "
                                  "invalid: regex mismatch", "ERROR");
      assert(m.size() == 5);

      assert(m[1].length());
      int minSlave = stoi(m[1]);
      int maxSlave = m[2].length() ? stoi(m[2]) : minSlave;
      if (maxSlave >= mpiSize)
        throw TCLAP::ArgException("given prefix MPI partition string is "
                                  "invalid: Slave " + to_string(maxSlave) +
                                  "exceeds MPI process count (" +
                                  to_string(mpiSize) + ")", "ERROR");
      if (minSlave == 0)
        throw TCLAP::ArgException("given prefix MPI partition string is "
                                  "invalid: Master 0 cannot be assigned to be "
                                  "a slave", "ERROR");
                   
      if (minSlave > maxSlave) 
        throw TCLAP::ArgException("given prefix MPI partition string is "
                                  "invalid: Slave numbering error: "
                                  "minSlave (" + to_string(minSlave) + ") "
                                  "exceeds maxSlave(" + to_string(maxSlave) +
                                  ")", "ERROR");
      
      assert(m[3].length());
      auto minLvl = std::stoul(m[3]);
      auto maxLvl = m[4].length() ? std::stoul(m[4]) : minLvl;

      if (minLvl >= args.prefix.size()) 
        throw TCLAP::ArgException("given prefix MPI partition string is "
                                  "invalid: Level numbering error: "
                                  "minLvl (" + to_string(minLvl) + ")"
                                  " exceeds prefix size (" +
                                  to_string(args.prefix.size()) + ")", "ERROR");

      if (maxLvl >= args.prefix.size()) 
        throw TCLAP::ArgException("given prefix MPI partition string is "
                                  "invalid: Level numbering error: "
                                  "maxLvl (" + to_string(maxLvl) + ")"
                                  " exceeds prefix size (" +
                                  to_string(args.prefix.size()) + ")", "ERROR");

      if (minLvl > maxLvl) 
        throw TCLAP::ArgException("given prefix MPI partition string is "
                                  "invalid: Level numbering error: "
                                  "minLvl (" + to_string(minLvl) + ") "
                                  "exceeds maxLvl (" + to_string(maxLvl) + ")",
                                  "ERROR");

      for (auto i = minLvl; i <= maxLvl; ++i)
        for (auto j = minSlave; j <= maxSlave; ++j)
          slavesForLvl[i].push_back(j);
    }

    for (size_t lvl = 0; lvl < args.prefix.size(); ++lvl) {
      for (auto rank : slavesForLvl[lvl])
        lvlsForSlave[rank].push_back(lvl);
    }

    for (size_t lvl = 0; lvl < args.prefix.size(); ++lvl) {
      if (slavesForLvl[lvl].empty())
        throw TCLAP::ArgException("Level " + to_string(lvl) + " has no "
                                  "slaves!");
    }

    for (int rank = 1; rank < mpiSize; ++rank) {
      if (lvlsForSlave[rank].empty())
        throw TCLAP::ArgException("Slave " + to_string(rank) + " has no "
                                  "levels!");
    }

    if (!lvlsForSlave[0].empty())
      throw TCLAP::ArgException("Master 0 has levels!");

    if (mpiRank > 0)
      args.workStack = new reduce::HierarchicalSlaveStack(lvlsForSlave[mpiRank],
                                                          args.prefix.size(),
                                                          mpiBufferSize);
    else
      args.workStack = new reduce::HierarchicalMasterStack(slavesForLvl);
  }
  else if (mpiRank == 0) {
#endif // WITH_OPEN_MPI
    args.workStack = new reduce::SoloStack();
#ifdef WITH_OPEN_MPI
  }
  else if (mpiRank > 0) {
    args.workStack = new reduce::VerySimpleSlaveStack(mpiBufferSize);
  }
  else {
    assert (false && "Not implemented yet!");
  }
#endif // WITH_OPEN_MPI

  if (estimationValue.getValue() != "") {
    vector<double> probs;
    std::istringstream iss(estimationValue.getValue());
    string s;
    while (std::getline(iss, s, ',')) {
      double d = std::stod(s);
      if (d < 0 || d > 1)
        throw TCLAP::ArgException("Token `" + s + "' violates range "
                                  "constraint [0,1]!");
      args.estimationProbabilities.push_back(d);
    }
    if (args.estimationProbabilities.size() != args.prefix.size())
      throw TCLAP::ArgException("Invalid number of probabilities!");
  }

  args.threshold = thresholdValue.getValue();
  args.symmetryOnly = symmetryOnlySwitch.getValue();
  args.incremental = incrementalSwitch.getValue();
  args.verbose = verboseSwitch.getValue();
  args.printAutSizes = printAutSizesSwitch.getValue();
  args.checkConflicts = conflictSwitch.getValue();
  args.logFilename = logFileValue.getValue();
#ifdef WITH_OPEN_MPI
  size_t p;
  while ((p = args.logFilename.find("%r")) != string::npos) {
    args.logFilename.replace(p, 2, to_string(mpiRank));
  }
#endif // WITH_OPEN_MPI
  
  return args;
}

static void output(const Args& args, const reduce::Reducer& r,
                   const vector<reduce::Reducer::Assignment>& assignments) {
  if (!args.symmetryOnly) {
    if(!args.incremental) {
      if(args.inCnf.empty()) {
        int count = 0;
        for (const auto& a : assignments) {
          assert(a.vars.size() == a.vals.size());
          count++;
          cout << count << ": [" << a.autSize << "] ";
          
          r.printAssignment(cout, a);
        }
      }
      else {
        /* Store conjuncts in a buffer. */
        vector<int> conjbuf;
        int cursor = 0;
        int count = 0;
        for (const auto& a : assignments) {
          assert(a.vars.size() == a.vals.size());
          count++;
          int len = a.vars.size(); 
          cerr << "c branch " << count << " " << a.autSize << endl;
          for(int i = 0; i < len; i++) {
            conjbuf.push_back((// a[1+i+len] ==
                               a.vals[i] ==
                               static_cast<int>(r.getVarMap().getValue(0))) ?
                              -(1+r.getVarMap().graphToCnf(a.vars[i])) :
                              1+r.getVarMap().graphToCnf(a.vars[i]));
            cursor++;
          }
          conjbuf.push_back(0);
          cursor++;
        }
        
        /* Print CNF with adjust for conjunct-clauses. */
        if (count == 0) {
          r.printCnf(cout, 
                      "cnf", 
                      count, 
                      0);
        }
        else {
          // ASSUMING THERE IS ANY
        
          r.printCnf(cout, 
                      "cnf", 
                      count, 
                      cursor - count + 1);
          /* Print the conjunct-clauses. */
          int nv_base = r.getClauses().getVariableCount();
          int u = 0;
          int end = cursor;
          cursor = 0;
          while(cursor < end) {
            if(conjbuf[cursor] == 0) {
              u++;
            }
            else {
              cout << conjbuf[cursor] << " " << -(1 + nv_base + u)
                   << " 0" << endl;
            }
            cursor++;
          }
          if(u != count)
            throw std::runtime_error("bad conjunct buffer");
          /* Print the final clause of conjunct-variables. */
          for(int i = 0; i < count; i++) {
            cout << (1 + nv_base + i) << (i == count - 1 ? " 0\n" : " ");
          }
        }
      }
    }
    else {
      r.printCnf(cout, "inccnf", -1, -1);
      int count = 0;
      for (const auto& a : assignments) {
        assert(a.vars.size() == a.vals.size());
        count++;
        cerr << "c branch " << count << " " << /* a[2*a[0]+1] */ a.autSize << endl;
        for(size_t i = 0; i < a.vars.size(); i++) {
          cout << (i == 0 ? "a " : " ")
               << (a.vals[i] == static_cast<int>(r.getVarMap().getValue(0)) ?
             -(1+r.getVarMap().graphToCnf(a.vars[i])) :
                   1+r.getVarMap().graphToCnf(a.vars[i]));
        }
        cout << " 0" << endl;
      }
    }
  }

  enable_timing(); // enable timings
  cerr << "host: " << common_hostname();
  pop_print_time("total");
  cerr << endl;

  cerr << "build: " << COMMITID << endl;
}

int main(int argc, char** argv) {
#ifdef WITH_OPEN_MPI
  MPI_Init(&argc,&argv);

  int mpiRank, mpiSize, mpiVersionLen;
  char mpiVersion[MPI_MAX_LIBRARY_VERSION_STRING];

  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Get_library_version(mpiVersion, &mpiVersionLen);
  
  cerr << "MPI init " << (mpiRank ? "slave " : "master ")
       << (mpiRank + 1) << " / " << mpiSize
       << ", version " << mpiVersion <<  endl;
#endif // WITH_OPEN_MPI
  
  Args args;
  try {
    args = parseArgs(argc,argv);
  }
  catch (TCLAP::ArgException& e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  enable_timing(); // enable timings

  push_time();
  push_time(); 

  reduce::Reducer r(args.verbose, args.printAutSizes,
                    args.checkConflicts,
                    args.inCnf, args.prefix,
                    args.threshold, args.graph,
                    args.varMap, args.workStack,
                    args.estimationProbabilities
#ifdef WITH_OPEN_MPI
                    , args.lvlsForSlave
#endif // WITH_OPEN_MPI
		    );

  disable_timing(); // time only the init phase

  auto beforeRun = std::chrono::steady_clock::now();
  
  vector<reduce::Reducer::Assignment> assignments;
  r.getPrefixAssignments(assignments);

  auto afterRun = std::chrono::steady_clock::now();
  std::chrono::nanoseconds diff = afterRun - beforeRun;

  std::ofstream ofs;
  if (args.logFilename != "")
    ofs.open(args.logFilename);
  std::ostream& logStream = (args.logFilename == "" ? std::cerr : ofs);
 
#ifdef WITH_OPEN_MPI
  if (mpiRank == 0)  {
#endif // WITH_OPEN_MPI
    output(args, r, assignments);  
#ifdef WITH_OPEN_MPI
  }  
  reduce::mpi::printMpiStatistics(logStream);
#endif // WITH_OPEN_MPI

  logStream << "total runtime: " << diff.count() << " ns" << endl;
  logStream << "conflict count: " << r.getConflictCount() << endl;
  
  logStream << "total nauty time: " << reduce::getTotalNautyTime() << " ns"
            << endl;;
  logStream << "number of nauty calls: " << reduce::getNautyCallCount() << endl;
  char hostname[1024];
  gethostname(hostname, 1024);
  logStream << "hostname: " << hostname << endl;
  char temp[256];
  sprintf(temp, 
          "c %7s %14s %14s %14s\n",
          "Size",
          "Generated",
          "Canonical",
          "Output");
  logStream << temp;
  
  for (size_t l = 0; l < r.getPrefixLength(); ++l)
    r.printStats(logStream, l);

#ifdef WITH_OPEN_MPI
  MPI_Finalize();
#endif // WITH_OPEN_MPI

  return EXIT_SUCCESS;
}
