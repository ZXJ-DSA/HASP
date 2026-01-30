
#include "head.h"

#include <boost/program_options.hpp>

#define DEFAULT_INDEX DVPLIndex
#define DEFAULT_PNUM 32
#define DEFAULT_RNUM 4
#define DEFAULT_QUERY 4
#define DEFAULT_PMETHOD "NC"
#define DEFAULT_UTYPE 0
#define DEFAULT_BSIZE 100
#define DEFAULT_BNUM 1
#define DEFAULT_UINTERVAL 300
#define DEFAULT_RESPONSETIME 1.0
#define DEFAULT_QNUM 10000
#define DEFAULT_THREADNUM 15
#define DEFAULT_WORKNUM 1
#define DEFAULT_CACHESIZE 30
#define DEFAULT_BANDWIDTH 100
#define DEFAULT_LOWERB 0.1
#define DEFAULT_UPPERB 2
#define DEFAULT_PRETASK 0
namespace po = boost::program_options;

HTSPIndex validateAndConvertHTSPIndex(int inputIndex);
string getHTSPIndexName(HTSPIndex index);
bool validateParameters(const string &algoParti, int algoQuery, int updateType, int preTask);

int main(int argc, char **argv) {
    string DesFile = "./data/";
    string dataset = "NY";
    HTSPIndex systemIndex = HTSPIndex::DVPL;
    int algoChoice = static_cast<int>(systemIndex);
    int partitionNum = DEFAULT_PNUM;
    int algoQuery = DEFAULT_QUERY;
    string algoParti = DEFAULT_PMETHOD;
    int updateType = DEFAULT_UTYPE;
    int runtimes = DEFAULT_QNUM;
    int batchNum = DEFAULT_BNUM;
    int batchSize = DEFAULT_BSIZE;
    int batchInterval = DEFAULT_UINTERVAL;
    int threadNum = DEFAULT_THREADNUM;
    int preTask = DEFAULT_PRETASK;
    int regionNum = DEFAULT_RNUM;
    int bandwidth = DEFAULT_BANDWIDTH;
    double lowerB = DEFAULT_LOWERB;
    double upperB = DEFAULT_UPPERB;
    double T_r = DEFAULT_RESPONSETIME;//average query response time, s
    int workerNum = DEFAULT_WORKNUM;
    int cacheListSize = DEFAULT_CACHESIZE;

    try {
        // ===================== 1. Define Parameter Options (Original Logic) =====================
        po::options_description desc("Program Usage");
        desc.add_options()
                // Required parameters (no default value)
                ("source_path,s", po::value<string>()->required(),
                 "Source path (required), e.g. /export/project/xzhouby")
                ("dataset,d", po::value<string>()->required(), "Name of dataset (required), e.g. NY")
                ("index,i", po::value<int>()->required(),
                 "HTSP system index (required):\n"
                 "  3 = PMHL\n"
                 "  5 = PostMHL\n"
                 "  6 = DVPL\n"
                 "  7 = DVPLwoR\n"
                 "  8 = DVPLwoRV\n"
                 "  9 = DVPLwoRVO\n"
                 "  10 = MVCC_PostMHL")
                // Optional parameters (with original default values & variable names)
                ("response_time,r", po::value<double>()->default_value(DEFAULT_RESPONSETIME),
                 "Average query response time requirement (seconds), default: 1.0")
                ("partition_num,p", po::value<int>()->default_value(DEFAULT_PNUM), "Partition number, default: 32")
                ("partition_method,m", po::value<string>()->default_value(DEFAULT_PMETHOD),
                 "Partition method (NC: PUNCH; MT: METIS), default: NC")
                ("query_strategy,q", po::value<int>()->default_value(DEFAULT_QUERY),
                 "Query strategy (0:A*; 1: PCH; 2: No-boundary; 3: Post-boundary; 4: Cross-boundary), default: 4")
                ("update_type,u", po::value<int>()->default_value(DEFAULT_UTYPE),
                 "Update type (0: No Update Test; 1: Decrease; 2: Increase), default: 0")
                ("batch_num,b", po::value<int>()->default_value(DEFAULT_BNUM), "Batch number, default: 10")
                ("cache_list_size,c", po::value<int>()->default_value(DEFAULT_CACHESIZE), "Cache list number, default: 30")
                ("batch_interval,I", po::value<int>()->default_value(DEFAULT_UINTERVAL), "Batch interval (seconds), default: 300")
                ("thread_num,T", po::value<int>()->default_value(DEFAULT_THREADNUM), "Thread number, default: 15")
                ("worker_num,W", po::value<int>()->default_value(DEFAULT_WORKNUM), "Query worker number, default: 1")
                ("region_num,R", po::value<int>()->default_value(DEFAULT_RNUM), "Region number for VPL, default: 4")
                ("bandwidth,B", po::value<int>()->default_value(DEFAULT_BANDWIDTH), "Bandwidth, default: 100")
                ("lower_bound_ratio,L", po::value<double>()->default_value(DEFAULT_LOWERB), "Lower bound ratio, default: 0.1")
                ("upper_bound_ratio,U", po::value<double>()->default_value(DEFAULT_UPPERB), "Upper bound ratio, default: 2")
                ("preprocessing_task,P", po::value<int>()->default_value(DEFAULT_PRETASK),
                 "Preprocessing task (1: Tree height-aware PSP Ordering; 2: Partitioned Query Generation), default: 0")

                // Help information
                ("help,h", "Display this help message");

        // ===================== 2. Parse Command Line Parameters =====================
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        // Display help message
        if (vm.count("help")) {
            cout << desc << endl;
            return 0;
        }

        // Validate required parameters (throw exception if missing)
        po::notify(vm);

        // ===================== 3. Read Parameters (100% Original Variable Names) =====================
        // Required parameters (original variable names)
        DesFile = vm["source_path"].as<string>();       // Original: DesFile (source path)
        dataset = vm["dataset"].as<string>();           // Original: dataset
        algoChoice = vm["index"].as<int>();           // Original: algoChoice (HTSP index)
        systemIndex = validateAndConvertHTSPIndex(algoChoice);
        // Optional parameters (strictly original variable names)
        T_r = vm["response_time"].as<double>();         // Original: T_r (response time)
        partitionNum = vm["partition_num"].as<int>();      // Original: partitionNum
        algoParti = vm["partition_method"].as<string>();// Original: algoParti
        algoQuery = vm["query_strategy"].as<int>();        // Original: algoQuery
        updateType = vm["update_type"].as<int>();          // Original: updateType
        batchNum = vm["batch_num"].as<int>();              // Original: batchNum
        cacheListSize = vm["cache_list_size"].as<int>();   // Original: cacheListSize
        batchInterval = vm["batch_interval"].as<int>();    // Original: batchInterval
        threadNum = vm["thread_num"].as<int>();            // Original: threadNum
        workerNum = vm["worker_num"].as<int>();            // Original: workerNum
        regionNum = vm["region_num"].as<int>();            // Original: regionNum
        bandwidth = vm["bandwidth"].as<int>();             // Original: bandwidth
        lowerB = vm["lower_bound_ratio"].as<double>();
        upperB = vm["upper_bound_ratio"].as<double>();
        preTask = vm["preprocessing_task"].as<int>();

        // ===================== 4. Validate Parameter Legality =====================
        if (!validateParameters(algoParti, algoQuery, updateType, preTask)) {
            return 1;
        }

        // ===================== 5. Output (Original Format + Original Variables) =====================
        cout << "argv[1] (Source Path): " << DesFile << endl;
        cout << "argv[2] (Dataset): " << dataset << endl;
        cout << "argv[3] (HTSP System Index): " << algoChoice << endl;
        cout << "argv[4] (Query Response Time, s): " << T_r << endl;
        cout << "argv[5] (Partition Number): " << partitionNum << endl;
        cout << "argv[6] (Partition Method): " << algoParti << endl;
        cout << "argv[7] (Query Strategy): " << algoQuery << endl;
        cout << "argv[8] (Update Type): " << updateType << endl;
        cout << "argv[9] (Batch Number): " << batchNum << endl;
        cout << "argv[10] (Cache List Size): " << cacheListSize << endl;
        cout << "argv[11] (Batch Interval): " << batchInterval << endl;
        cout << "argv[12] (Thread Number): " << threadNum << endl;
        cout << "argv[13] (Query Worker Number): " << workerNum << endl;
        cout << "argv[14] (Region Number): " << regionNum << endl;
        cout << "argv[15] (Bandwidth): " << bandwidth << endl;
        cout << "argv[16] (Lower Bound Ratio): " << lowerB << endl;
        cout << "argv[17] (Upper Bound Ratio): " << upperB << endl;
        cout << "argv[18] (Preprocessing Task): " << preTask << endl;
    } catch (const po::required_option &e) {
        // Missing required parameters
        cerr << "Error: Missing required parameter!" << endl;
        cerr << e.what() << endl;
        return 1;
    } catch (const po::invalid_option_value &e) {
        // Invalid parameter value type (e.g., string for integer parameter)
        cerr << "Error: Invalid parameter value type!" << endl;
        cerr << e.what() << endl;
        return 1;
    } catch (const exception &e) {
        // Other exceptions
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    //used for running time calculation
    Timer tt0;
    tt0.start();
    string sourcePath = DesFile + "/" + dataset;
    string ODfile = sourcePath + "/" + dataset + ".query";
    string updateFile = sourcePath + "/" + dataset + ".update";

    Graph g;
    g.threadnum = threadNum;//thread number of parallel computation (can be changed)
    g.sourcePath = sourcePath;
    g.ifParallel = true;
    g.dataset = dataset;
    g.algoQuery = algoQuery;
    g.algoChoice = algoChoice;
    g.algoUpdate = algoQuery;
    g.algoParti = algoParti;
    g.partiNum = partitionNum;
    g.bandWidth = bandwidth;
    g.bRatioLower = lowerB;
    g.bRatioUpper = upperB;
//    g.samePartiPortion=portion;
    g.lambdaCacheSize = cacheListSize;
    cout << "Dataset: " << dataset << endl;
    cout << "System Index: " << getHTSPIndexName(systemIndex) << endl;
    if (g.algoQuery == Dijk) {
        cout << "Dijkstra's test !!!!!!!" << endl;
    } else if (g.algoQuery == PH2H_No) {
        cout << "This is test for No-boundary strategy !!!!!!!" << endl;
    } else if (g.algoQuery == PH2H_Post) {
        cout << "This is test for Post-boundary strategy !!!!!!!" << endl;
    } else if (g.algoQuery == PH2H_Cross) {
        cout << "This is test for Cross-boundary strategy !!!!!!!" << endl;
    } else if (g.algoQuery == PCH_No) {
        if (g.algoChoice == 2) {
            cout << "Wrong index choice." << endl;
            exit(1);
        }
        cout << "This is test for PCH-No !!!!!!!" << endl;
    } else {
        cout << "Wrong query strategy! " << g.algoQuery << endl;
        exit(1);
    }
    cout << "Partition method: " << g.algoParti << endl;
    cout << "Partition number: " << g.partiNum << endl;
    cout << "Thread number: " << g.threadnum << endl;
    if (workerNum > threadNum) {
        workerNum = threadNum;
    }
    cout << "Query worker number: " << workerNum << endl;

    if (preTask == 1) {
        g.PH2HVertexOrdering(2);//Naive Boundary-first MDE ordering
//        g.PH2HVertexOrdering(3);//Tree height-aware PSP ordering
    } else if (preTask == 2) {
        g.QueryGenerationParti(true);//same partition and real-world simulation
    }

//    g.ReadGraph(graphfile);//
//    g.StainingMethod(0);
//    g.getThroughputEvolveData("/Users/zhouxj/Documents/1-Research/1-Papers/0-My_Papers/5-HASP/result/log_exp5_1");
    ///Task 1: Index construction
    g.IndexConstruction();
//    g.PH2HIndexConstruct();
//    g.WriteCTIndex(graphfile);

    ///Task 2: Query processing
    g.CorrectnessCheck(100);
    g.EffiCheck(ODfile, runtimes);//query efficiency test
//    g.EffiCheck(ODfile+"Parti",runtimes);//query efficiency test
//    g.EffiCheck(ODfile+"SameParti",runtimes);
//    g.EffiCheck(ODfile+"CrossParti",runtimes);
//    g.MHLIndexCompareCH(g.sourcePath+dataset+".CHIndex0");
//    exit(0);
    ///Task 3: Index update
    g.RealUpdateThroughputTestQueueModel(
            sourcePath + "/" + dataset + "_20160105_" + to_string(batchInterval) + ".batchUpdates", batchNum, T_r,
            workerNum, regionNum);

    tt0.stop();
    cout << "\nOverall runtime: " << tt0.GetRuntime() << " s." << endl;
    cout << "------------------\n" << endl;
//    exit(0);
    g.clear();

    return 0;
}

// ===================== Core Parameter Validation Function (Original Variable Names) =====================
bool validateParameters(const string &algoParti, int algoQuery, int updateType, int preTask) {
    bool isValid = true;

    // Validate partition method (NC/MT)
    if (algoParti != "NC" && algoParti != "MT") {
        cerr << "Error: Partition method must be NC (PUNCH) or MT (METIS)! Input value: " << algoParti << endl;
        isValid = false;
    }

    // Validate query strategy (0-4)
    if (algoQuery < 0 || algoQuery > 4) {
        cerr << "Error: Query strategy must be 0-4! Input value: " << algoQuery << endl;
        isValid = false;
    }

    // Validate update type (0-2)
    if (updateType < 0 || updateType > 2) {
        cerr << "Error: Update type must be 0-2! Input value: " << updateType << endl;
        isValid = false;
    }

    // Validate preprocessing task (0-2)
    if (preTask != 0 && preTask != 1 && preTask != 2) {
        cerr << "Error: Preprocessing task must be 0, 1 or 2! Input value: " << preTask << endl;
        isValid = false;
    }

    return isValid;
}

HTSPIndex validateAndConvertHTSPIndex(int inputIndex) {
    // Check if input is in the range of legal enum values
    switch (inputIndex) {
        case static_cast<int>(HTSPIndex::PMHL):
            return HTSPIndex::PMHL;
        case static_cast<int>(HTSPIndex::PostMHL):
            return HTSPIndex::PostMHL;
        case static_cast<int>(HTSPIndex::DVPL):
            return HTSPIndex::DVPL;
        case static_cast<int>(HTSPIndex::DVPLwoR):
            return HTSPIndex::DVPLwoR;
        case static_cast<int>(HTSPIndex::DVPLwoRV):
            return HTSPIndex::DVPLwoRV;
        case static_cast<int>(HTSPIndex::DVPLwoRVO):
            return HTSPIndex::DVPLwoRVO;
        case static_cast<int>(HTSPIndex::MVCC_PostMHL):
            return HTSPIndex::MVCC_PostMHL;
        default:
            // Throw clear error (includes all legal values for user reference)
            throw invalid_argument(
                    "Invalid HTSP system index: " + to_string(inputIndex) +
                    "\nLegal values are:\n"
                    "  3 = PMHL\n"
                    "  5 = PostMHL\n"
                    "  6 = DVPL\n"
                    "  7 = DVPLwoR\n"
                    "  8 = DVPLwoRV\n"
                    "  9 = DVPLwoRVO\n"
                    "  10 = MVCC_PostMHL"
            );
    }
}

string getHTSPIndexName(HTSPIndex index) {
    switch (index) {
        case HTSPIndex::PostMHL: return "PostMHL";
        case HTSPIndex::DVPL: return "DVPL";
        case HTSPIndex::DVPLwoR: return "DVPLwoR";
        case HTSPIndex::DVPLwoRV: return "DVPLwoRV";
        case HTSPIndex::DVPLwoRVO: return "DVPLwoRVO";
        default: return "UNKNOWN"; // Should never reach here (due to validation)
    }
}