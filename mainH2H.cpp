#include "headH2H.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    // ===================== 1. Define Original Variables (Keep All Names) =====================
    string DesFile = "./data/";
    string dataset = "NY";
    int algorithm = 2;
    int updateType = 0;
    int runtimes = 1000; // Fixed value, no need to parse
    // runtimes=100; // Keep original comment
    int batchNum = 10;
    int batchSize = 10; // Original unused variable, keep
    int batchInterval = 60; // Original default value (not 10, follow variable init)
    int threadNum = 15;
    double T_r = 1; // Average query response time (default 1)
    string queryFName;
    int workerNum = 15;
    int cacheSize = 20;

    try {
        // ===================== 2. Define Boost Parameter Options =====================
        po::options_description desc("Program Usage");
        desc.add_options()
                // Required parameters (match original argv[1]-argv[3], no default)
                ("source_path,s", po::value<string>()->required(),
                 "Source path (required), e.g. /export/project/xzhouby")
                ("dataset,d", po::value<string>()->required(), "Name of dataset (required), e.g. NY")
                ("algorithm,a", po::value<int>()->required(),
                 "Algorithm (required): 0=Dijkstra; 1=CH; 2=H2H; 3=BiDijkstra")

                // Optional parameters (match original argv[4]-argv[11], keep original defaults)
                ("response_time,r", po::value<double>()->default_value(1),
                 "Average query response time requirement (seconds), default: 1")
                ("update_type,u", po::value<int>()->default_value(0),
                 "Update type (0: No Update Test; 1: Decrease; 2: Increase), default: 0")
                ("batch_num,b", po::value<int>()->default_value(10), "Batch number, default: 10")
                ("cache_list_size,c", po::value<int>()->default_value(20), "Cache list number, default: 20")
                ("batch_interval,i", po::value<int>()->default_value(60),
                 "Batch interval (seconds), default: 60") // Follow original variable default (60)
                ("thread_num,t", po::value<int>()->default_value(15), "Thread number, default: 15")
                ("worker_num,w", po::value<int>()->default_value(15), "Query worker number, default: 15")
                ("query_file,f", po::value<string>(), "Query file name (optional, no default)")

                // Help information
                ("help,h", "Display this help message");

        // ===================== 3. Parse Command Line Parameters =====================
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        // Display help message
        if (vm.count("help")) {
            cout << desc << endl;
            return 0;
        }

        // Validate required parameters (throw exception if missing)
        po::notify(vm);

        // ===================== 4. Read Parameters (100% Original Variable Names) =====================
        // Required parameters (match Boost option names → original variables)
        DesFile = vm["source_path"].as<string>();       // Original: argv[1] → DesFile
        dataset = vm["dataset"].as<string>();           // Original: argv[2] → dataset
        algorithm = vm["algorithm"].as<int>();          // Original: argv[3] → algorithm

        // Optional parameters (Boost auto uses default if not specified)
        T_r = vm["response_time"].as<double>();         // Original: argv[4] → T_r
        updateType = vm["update_type"].as<int>();       // Original: argv[5] → updateType
        batchNum = vm["batch_num"].as<int>();           // Original: argv[6] → batchNum
        cacheSize = vm["cache_list_size"].as<int>();    // Original: argv[7] → cacheSize
        batchInterval = vm["batch_interval"].as<int>(); // Original: argv[8] → batchInterval
        threadNum = vm["thread_num"].as<int>();         // Original: argv[9] → threadNum
        workerNum = vm["worker_num"].as<int>();         // Original: argv[10] → workerNum

        // Optional: query_file (no default, only set if specified)
        if (vm.count("query_file")) {
            queryFName = vm["query_file"].as<string>(); // Original: argv[11] → queryFName
        }

        // ===================== 5. Output Parameters (Original Format) =====================
        cout << "argc: " << argc << endl; // Keep original output
        cout << "argv[1] (Source Path): " << DesFile << endl;
        cout << "argv[2] (Dataset): " << dataset << endl;
        cout << "argv[3] (Algorithm): " << algorithm << endl;
        cout << "argv[4] (Query Response Time): " << T_r << endl;
        cout << "argv[5] (Update Type): " << updateType << endl;
        cout << "argv[6] (Batch Number): " << batchNum << endl;
        cout << "argv[7] (Cache List Number): " << cacheSize << endl;
        cout << "argv[8] (Batch Interval): " << batchInterval << endl;
        cout << "argv[9] (Thread Number): " << threadNum << endl;
        cout << "argv[10] (Query Worker Number): " << workerNum << endl;
        if (!queryFName.empty()) {
            cout << "argv[11] (Same-parti Proportion): " << queryFName << endl; // Keep original comment/label
        }

    } catch (const po::required_option &e) {
        // Missing required parameters (original argv[1]-argv[3])
        cerr << "Error: Missing required parameter!" << endl;
        cerr << e.what() << endl;
        return 1;
    } catch (const po::invalid_option_value &e) {
        // Invalid value type (e.g., string for integer parameter)
        cerr << "Error: Invalid parameter value type!" << endl;
        cerr << e.what() << endl;
        return 1;
    } catch (const exception &e) {
        // Other exceptions (e.g., bad cast, unknown option)
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    //used for running time calculation
    Timer tt0;
    tt0.start();

    string sourcePath = DesFile + "/" + dataset + "/";
    string ODfile = sourcePath + dataset + ".query";
    string updateFile = sourcePath + dataset + ".update";

    Graph g;
    g.dataset = dataset;
    g.threadnum = threadNum;
    g.sourcePath = sourcePath;
    g.algoIndex = algorithm;
    g.queryFName = queryFName;
    g.lambdaCacheSize = cacheSize;
    cout << "Dataset: " << dataset << endl;
    cout << "Thread number: " << threadNum << endl;
    if (workerNum > threadNum) {
        workerNum = threadNum;
    }
    cout << "Query worker number: " << workerNum << endl;


    g.ReadGraph(sourcePath + dataset);//
//    g.StainingMethod(0);

    ///Task 1: Index construction
    g.IndexConstruction(algorithm);
//    g.H2HIndexConstruct();
//    g.WriteCTIndex(graphfile);

    ///Task 2: Query processing
//    g.CorrectnessCheckCore(100);
    g.CorrectnessCheck(100);
    g.EffiCheck(runtimes);

//    g.SameTreeQueryTest(ODfile,runtimes);
//    exit(0);

    ///Task 3: Index update
//    if(dataset=="beijing" || dataset=="Guangdong"){//real-life updates
    g.RealUpdateThroughputTestQueueModel(
            sourcePath + dataset + "_20160105_" + to_string(batchInterval) + ".batchUpdates", batchNum, T_r, workerNum);
//    }else{
//        g.RandomUpdateThroughputTestQueueModel(batchNum, batchSize, batchInterval, T_r, workerNum);
//        g.RandomUpdateThroughputTest(sourcePath+dataset+".update", batchNum, batchSize, batchInterval);
//        g.SPThroughputTest(updateType, ifBatch, batchNum, batchSize, batchInterval, runtimes);
//    g.IndexMaintenanceCHWP(updateType, updateSize, ifBatch, batchSize);//index maintenance
//    g.IndexMaintenanceH2H(updateType, updateSize, ifBatch, batchSize);//index maintenance
//    g.IndexMaintenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance
//    }


    tt0.stop();
    cout << "\nOverall runtime: " << tt0.GetRuntime() << " s." << endl;
    cout << "------------------\n" << endl;
//    exit(0);
    g.clear();
    return 0;
}
