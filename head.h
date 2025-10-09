
#ifndef HEAD_H_
#define HEAD_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <boost/thread/thread.hpp>
#include <chrono>
#include <string>
#include <algorithm>
#include "RTree.h"
#include "Heap.h"
//#include "labeling.hpp"
//#include <omp.h>

//typedef RTree<Rect*, int, 2, int> MyTree;

#define PI 3.1415926
#define INF 99999999
#define Dijk 0
#define PCH_No 1
#define PH2H_No 2
//#define PH2H_overlay 2
#define PH2H_Post 3
#define PH2H_Cross 4
#define CH 1
#define H2H 4
#define MICROSEC_PER_SEC 1000000
#define MILLISEC_PER_SEC 1000
#define Resolution 10000000
#define DEC -1
#define INC 1
#define AFF 0
#define UNAFFIN 1
#define UNAFF 2
#define DecUpdate 1
#define IncUpdate 2
#define MixUpdate 3
//typedef unsigned int vertex;
typedef int vertex;

using namespace std;
//using namespace boost;


extern vector<int> NodeOrder_;//nodeID order
extern vector<int> _DD_;
extern vector<int> _DD2_;

struct Nei{
	int nid;
	int w;
	int c;
};

struct OrderCompMax{// maximum-first, Higher-order first
    int ID;
    OrderCompMax(){ID=0;}
    OrderCompMax(int _ID){
        ID=_ID;
    }
    bool operator< (const OrderCompMax d) const{
        if(NodeOrder_[ID]!=NodeOrder_[d.ID])
            return NodeOrder_[ID]>NodeOrder_[d.ID];
        return ID>d.ID;
    }
};

struct OrderCompMin{//prior to reture the vertex with smaller order
    int x;
    OrderCompMin(int _x){
        x=_x;
    }
    bool operator< (const OrderCompMin& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
        }
    }
};

struct OrderCompp{//prior to return the vertex with smaller order
    int x;
    OrderCompp(int _x){
        x=_x;
    }
    bool operator< (const OrderCompp& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
        }
    }
};


//return smallest Dis with largest-order vertex
struct MinComp{
    int s, t, Dis;//s is the source vertex while t is the target vertex
    MinComp(){
        s=0, t=0, Dis=0;
    }
    MinComp(int _ID1, int _ID2, int _Dis){
        s=_ID1; t=_ID2; Dis=_Dis;
    }
    bool operator< (const MinComp d) const{
        if(Dis != d.Dis){
            return Dis<d.Dis;
        }else{
            if(NodeOrder_[s]!=NodeOrder_[d.s]){
                return NodeOrder_[s]>NodeOrder_[d.s];
            }
            return NodeOrder_[t]>NodeOrder_[d.t];
        }
    }
};

//return the largest value
struct MaxComp{// maximum-first,
    int value;
    int pid1=-1,pid2=-1;
    MaxComp(){value=0;}
    MaxComp(int _value, int _pid1, int _pid2){
        value=_value, pid1=_pid1, pid2=_pid2;
    }
    bool operator< (const MaxComp d) const{
        if(value!=d.value){
            return value>d.value;
        }else{
            if(pid1!=d.pid1){
                return pid1>d.pid1;
            }
            return pid2>d.pid2;
        }

    }
};

struct DegComp{//min-first
    int x;
    DegComp(int _x){
        x=_x;
    }
    bool operator < (const DegComp d) const{
        if(_DD_[x]!=_DD_[d.x])
            return _DD_[x]<_DD_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct DegHeightComp{//min-first
    int x;
    DegHeightComp(int _x){
        x=_x;
    }
    bool operator < (const DegHeightComp d) const{
        if(_DD_[x]!=_DD_[d.x]){
            return _DD_[x]<_DD_[d.x];
        }else{
            if(_DD2_[x]!=_DD2_[d.x]){
                return _DD2_[x]<_DD2_[d.x];
            }else{
//                return true;
                return x<d.x;
            }
        }
    }
};

struct DegComp2{//min-first
    int x;
    DegComp2(int _x){
        x=_x;
    }
    bool operator < (const DegComp2 d) const{
        if(_DD2_[x]!=_DD2_[d.x])
            return _DD2_[x]<_DD2_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count(how many ways can lead to this super edge weight)
	vector<int> pos;
	vector<int> dis, cnt;//the distance value and corresponding count number (i.e., record how many path has the shortest distance)
//    unordered_map<int,int> disOldM;//record the old distance for update
    vector<int> disPost, cntPost;//the distance value of post-boundary strategy and corresponding count number
    vector<int> disInf;//the distance from this vertex to the boundary vertices and corresponding count number, cntInf
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
    int rootPos;//index position of the root vertex of this partition in vAncestor
	vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert), crucial for prune invalid neighbor check
    vector<bool> FNPost;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert)
    vector<bool> FNInf;//whether the interface distance is obtained from shortcuts (vert)
    vector<int> disChange;//indicate whether the distance to the ancestors have changed, >0: increase, <0: decrease, 0: no update
    bool ifU=false;//indicate whether need to update
//    vector<bool> vertDisRe;//indicate whether the shortcuts to the vert neighbors have changed
	unordered_set<int> DisRe;//record the vert vertex id that the distance label should be updated
    unordered_set<int> DisRePost;//record the vert vertex id that the interface label should be updated
//    unordered_set<int> disDec;//record the ancestor whose label should be decreased
	vector<int> ch;//position index of children
	int height, hdepth;//hdepty is the deepest node that a vertex still exists
	int pa;//parent, the pa of root vertex is 0, position index
	int uniqueVertex;//vertex id of this tree node
//	vector<int> piv;//pivot vetex, used in path retrieval
//    int treeroot;//the tree id of subtree root, i.e., rank[x]
	Node(){
		vert.clear();
		pos.clear();
		dis.clear(); cnt.clear();
        disPost.clear(); cntPost.clear();
        disInf.clear();
        vAncestor.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
		FN.clear(); FNInf.clear();
		DisRe.clear();
        DisRePost.clear();
//		piv.clear();
//        treeroot=-1;
	}
};


struct simulate_times{
    std::vector<std::vector<int> > query_times;
    std::vector<int> update_times;
    int mode=2; //1: query_first 2: update_first
};

/// Throughput related function and data structure
struct query{
    double init_time=0; // in seconds, time of arrival
    double process_time=0;  // in seconds, query processing time
    query(double initT, double processT) {
        init_time = initT;
        process_time = processT;//query processing time
    }
};


class Graph{
public:
    string sourcePath;// the source path
    string dataset;
	int node_num=0;    //vertex number
	unsigned long long edge_num=0;    //edge number
	vector<vector<pair<vertex,int>>> Neighbor;//original graph
    vector<vector<pair<vertex,int>>> NeighborsParti;//<node_number,<in-partition adjacency lists>>
    vector<unordered_map<vertex,int>> NeighborsPartiPost;//adjacency lists of post-boundary partitions
    vector<unordered_map<vertex,int>> NeighborsOverlay;//<node_number,<adjacency lists of overlay graph>>
    vector<vector<pair<vertex,int>>> NeighborsOverlayV;//<node_number,<adjacency lists of overlay graph>>
//    vector<unordered_map<int,int>> BoundaryShortcuts;
    vector<pair<int,bool>> PartiTag;//<node_number,<partition_id,if_boundary>>, for PMHL
    vector<pair<int,set<int>>> PartiTags;//<node_number,<flag,set<partition_id>>>, flag=-1 means boundary vertex, flag>=0 means partition id
    vector<vector<vertex>> PartiVertex;//<partition_number,<in-partition vertices>>, in increasing vertex order, higher-rank vertex first
    vector<vector<vertex>> BoundVertex;//boundary vertices of each partition
//    vector<unordered_set<int>> BoundVertexSet;//boundary vertices of each partition
    vector<unordered_map<int,int>> BoundVertexMap;//map from boundary vertex id to its position id in BoundVertex
    vector<vertex> OverlayVertex;//overlay vertex in decreasing vertex order
//    vector<unordered_map<vertex,pair<int,int>>> repairShortcuts;//<ID1,ID2,<overlay weight, in-partition weight>>
    vector<map<int,int>> BoundShortcuts;//<ID1,ID2,weight>, source=-1 means stemming from core while source>=0 means stemming from specific partition
    int partiNum;   //partition number
    bool ifParallel = true;
    bool ifIncrease = false;
    int nodeNumOverlay = 0;//vertex number of overlay graph
    vector<double> uStageDurations;//durations of update stages
    vector<double> stageDurations;//durations of query stages
    vector<pair<int,int>> GraphLocation;//Longitude and Latitude, already times 1000,000

    //PostMHL
    double bRatioLower=0.2;//imbalance ratio of partitions
    double bRatioUpper=2;
    int bandWidth=100;//50;//bandwidth of PostMHL
    int HighestOrder;
    vector<bool> existOverlay;
    vector<int> partiRoots;//vertex id of partition root
    vector<vector<int>> partiRootsV;//vertex id of partition root
    vector<map<int,map<int,int>>> SuppPartiID;//<ID1,<ID2,<pid,weight>>>, record the partition and its supportive value for a given interface edge
    vector<map<int,pair<int,set<int>>>> SuppPartiIDReal;//ID1,<ID2,<weight,set<pid>>>>, //record the partitions that really support a given interface edge
    set<int> affectedParti;//record the affected partitions
    set<int> affectedPartiInc;//record the affected partitions
    vector<int> ProBeginIDV;
    vector<int> childNums;//children number of each vertex

	vector<int> DD; //intermediate variable in Contraction, DD2
	int threadnum=15;  //thread number
    int algoQuery=0;//algorithm for querying. 0:Dijkstra; 1: No-boundary; 2: Post-boundary; 3: Extended Label. 0: Dijkstra; 1: CH; 2:H2H. 0: Dijkstra; 1: PCH-No; 2: PH2H-No; 3: PCH-Post; 4: PH2H-Post; 5: PH2H-Extend
    int algoUpdate=0;//algorithm for core construction, (0: BPCL; 1: PCL; 2: PLL; 3: WPSL; 4: GLL; 5: Read)
    string algoParti="NC";
    int algoChoice=1;//HTSP system index. 1: Non-partition; 2: Partition

	//vertex order
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID
    int treewidth=0;//treewidth

    double tAncestor=0,tBoundary=0;

    /// Index Construction
//    vector<omp_lock_t> oml;
    unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
    Semaphore* sm = new Semaphore(1);;// = new Semaphore(threadnum);

    //H2H index construction
    //intermediate variable and function used in the H2H index construction
    vector<map<int,pair<int,int>>> E;//ID1,ID2,(weight,count)
//    vector<map<OrderCompMin, pair<int,int>>> EOrder;
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;//ID1,ID2,(weight,count)
    //for overlay graph
    vector<Node> Tree;
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//?
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
//    vector<int> EulerSeq;//prepare for the LCA calculation, EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    //    vector<map<int, vector<int>>> SCconNodesMT;//supportive vertex, multiple thread of SCconNodes
    vector<map<int, vector<pair<int,int>>>> SCconNodesMT;//<ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNid;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for overlay

    //for no-boundary partitions
    vector<vertex> IDMap;//map the old id to new id in partitions
    vector<vector<Node>> Trees;//Trees for no-boundary
    vector<vector<int>> toRMQs;
    vector<vector<vector<int>>> RMQIndexs;
    vector<vector<int>> ranks;
    vector<int> heightMaxs;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMTP;//for partitions. <ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNidP;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for partition
    int cliqueBoundaryPNum=0;
    //for post-boundary partitions
//    vector<vertex> IDMapPost;//map the old id to new id in partitions
    vector<vector<Node>> TreesPost;//Trees for post-boundary
    vector<vector<int>> toRMQsPost;
    vector<vector<vector<int>>> RMQIndexsPost;
    vector<vector<int>> ranksPost;
    vector<int> heightMaxsPost;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMTPost;//for partitions. <ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNidPost;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for partition

    vector<bool> ifRepaired;

    vector<int> ProBeginVertexSetOverlay;
    vector<int> ProBeginVertexSetOverlayInc;//for increase update of VPL
    set<int> vertexIDChLOverlay;
    vector<vector<int>> ProBeginVertexSetParti;
    vector<vector<int>> ProBeginVertexSetPartiInc;
    vector<set<int>> vertexIDChLParti;
    vector<vector<int>> ProBeginVertexSetPartiExtend;//for post-boundary update
    vector<vector<int>> ProBeginVertexSetPartiExtendInc;//for post-boundary update
//    vector<set<int>> vertexIDChLPartiExtend;

    //extension label
    vector<Node> TreeExt;
    vector<int> toRMQExt;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndexExt;//?
    vector<int> rankExt;//rank[v]>0 indicates non-core vertex
    int heightMaxExt;

    vector<bool> vUpdated;// flag of whether the label of vertex has been updated
    vector<int> bHeights;
    int samePartiPortion=-1;

    // for PSP ordering
    //    vector<pair<pair<int,int>,int>> CutEdges;//the cut edges
    vector<vector<int>> NeighborSketch;
    vector<unordered_set<int>> NeighborSketchS;
    vector<unordered_map<int,vector<int>>> NeighborSketchB;//<PID1,<PID2,<boundary IDs of PID2>>> the sketch graph and corresponding boundary vertex
//    vector<map<int,int>> vNodeOrderParti;
    vector<vector<int>> vNodeOrderParti;
    vector<int> vNodeOrderOverlay;
    double overlayUpdateT=0;

    //For VPL
    vector<unordered_map<int,unordered_set<int>>> BP_Flag;//<ID1, <ID2, PIDs>>
//    vector<map<pair<int,int>,unordered_set<int>>> BP_Flag;//for each partition, map from boundary edges (smaller ID, larger ID) to its affected partitions
    vector<vector<pair<int,int>>> PB_Flag;//for each partition, stores the boundary edges (smaller ID, larger ID) that lie in the shortest paths for its same-partition queries

    vector<bool> innerAffectedParti;//record whether a partition is inner-affected
    vector<bool> outerAffectedParti;//record whether a partition is outer-affected
    vector<bool> AffectedParti;//record whether a partition is affected
    int innerBenefitNum=0;
    int outerBenefitNum=0;
    int queryNum=0;//total query number before last query stage
    vector<pair<int,int>> overlayShortcutDec;//record the overlay shortcut updates caused by decrease update

    vector<Rect> partiMBR;//MBR of all partitions
    RTree<Rect*, int, 2, double> rtree;

    vector<pair<int,int>> pathComp1;
    vector<pair<int,int>> pathComp2;
    int updateBatchID=0;
    map<pair<int,int>,pair<int,int>> overlayLDec;//storing the source labels of decrease update
    map<pair<int,int>,int> overlayLInc;//storing the source labels of increase update
    vector<map<pair<int,int>,pair<int,int>>> partiLDec;//storing the source labels of decrease update
    vector<map<pair<int,int>,int>> partiLInc;//storing the source labels of increase update
    vector<long long int> checkNumParti;//record the check number for partition label update
    vector<vector<int>> clusterP;//record the regions, each contains several partitions
    vector<vector<double>> VPLUpdateTimes;//to record the detailed update time for VPL, <schedule,<U-Stages>>
    vector<vector<pair<bool,int>>> partiAffectInfo;//to record the affect info for each partition, <schedule,<parti,<V-Stage 1, V-stage 2(0: affected, 1: inner-unaffected, 2: unaffected)>>>
    unordered_map<int,vector<query>> lambdaCache; //query list cache for different lambda;
    vector<int> lambdaCount;
    benchmark::heap<2,int,int> lambdaPQueue;//lambda and its appear number
    int lambdaCacheSize=20;
    vector<double> unaffectedQNum;
    vector<double> totalQNum;


    ~Graph(){
        clear();
    }
    void clear(){
        Neighbor.clear();
        Tree.clear();
        vSm.clear();
    }
    //Dijkstra
    int Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    pair<int,vector<pair<int,int>>> DijkstraPath(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int BiDijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int Astar(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int EuclideanDis(int s, int t);
    void RetrievePath(int ID1, int ID2, vector<int> & prece, vector<pair<int,int>>& pathInfo);
    void RetrievePathCH(int ID1, int ID2, int finalID, vector<int> & preceF, vector<int> & preceB, vector<pair<int,int>>& pathInfo);
    int DijkstraCore(int ID1, int ID2);
    int DijkstraCorePath(int ID1, int ID2);


    /// For system throughput simulation
    pair<double,double> ThroughputEstimate(vector<vector<double>> &query_costs, vector<vector<double>> &update_costs, double threshold_time, double T);
    double analytical_update_first(double T_q, double T_u, double T_r, double T, double V_q);
    double ThroughputSimulate(vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime, int workerNum);
    double ThroughputSimulateReal(vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime, int workerNum);
    double ThroughputSimulateRealVPL(vector<pair<int,int>>& ODpair, vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime, int workerNum);
    pair<double, double> simulator_UpdateFirst(vector<vector<double>> & query_costs, vector<vector<double>> & update_costs, vector<query> &queryList, int T, double period_time);
    pair<double, double> simulator_UpdateFirst(vector<vector<double>> & query_costs, vector<vector<double>> & update_costs, vector<query> &queryList, int T, double period_time, int workerNum);
    pair<double, double> simulator_UpdateFirstVPL(vector<pair<int,int>>& ODpair, vector<vector<double>> & query_cost, vector<vector<double>> & update_cost, vector<query> &queryList, int T, double period_time, int workerNum);
    pair<double,int> getFastestWorker(vector<double>& workers);

    /// For non-partition SP index
    void HybridSPIndexConstruct();
    void MDEContraction(string orderfile);
//    void deleteEOrderGenerate(int u,int v);
    void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p);
    void insertEMTOrderGenerate(int u,int v,int w);
    //for contraction
    void deleteEorder(int u,int v);
    void insertEorder(int u,int v,int w);
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
    void makeTree();
    void makeIndex();
    void IndexsizeH2H();  //Index size computation

    /// For partitioned SP index
    void HybridPSPIndexConstruct();
    void PMHLIndexConstruct();
    void PMHLIndexConstructOpt();
    void PostMHLIndexConstruct();

    /// For PH2H
    void PH2HIndexConstruct(); //PH2H index construction
    void ConstructBoundaryShortcutV(vector<int> & p, bool ifAllPair);
    void ConstructBoundaryShortcut(int pid);
    void ConstructBoundaryShortcutNoAllPair(int pid);
    void Construct_PartiIndex(bool ifParallel, bool ifLabelC);
    void Construct_OverlayGraph(bool ifParallel);
    void Construct_OverlayGraphNoAllPair(bool ifParallel);
    void Construct_OverlayIndex(bool ifLabelC);

    void PreConstructAllPairs(bool ifParallel);
    void PreConstructAllPairsPartiV(vector<int> & p);
    void PreConstructAllPairsParti(int pid);

    /// Post-boundary index
    void RefreshBoundaryEdgesAndLabelingPartiV(vector<int>& p);
    void RefreshBoundaryEdgesAndLabelingParti(int pid);
    void ConstructPartitionPost(bool ifParallel);
    void ConstructPostParti(int pid);
    void ConstructPostPartiV(vector<int>& p);
    void ConstructPartitionPostIndex(bool ifParallel, bool ifLabelU);
    void ConstructPartitionPostIndexOpt(bool ifParallel);
    void Repair_PartiIndex(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void Repair_PartiIndexPostMHLPost(bool ifParallel, int updateType, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, double & runT);
    void Repair_PartiIndexForOpt(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, bool ifOpt);
    void RepairPartitionIndexDecrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexDecreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexIncrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexIncreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);

    void PostMHLIndexUpdatePostPartV(vector<int>& p, bool ifIncrease);
    void PostMHLIndexUpdatePostPart(int pid, bool ifIncrease);
    void PostMHLIndexUpdatePostPartDFS(int pid, int p, vector<int>& ancestor, vector<int> & interface, vector<map<int,int>>& disInfs, int rootHeight, bool ifIncrease);
    void PostMHLIndexUpdatePostPartDFS(int pid, int p, vector<int>& ancestor, vector<int> & interface, vector<map<int,int>>& disInfs, set<int>& vertexIDChL, int rootHeight, bool ifIncrease, Timer& tt);

    /// Extension index
    void makeTreeExtDFS(int ID, vector<int> &list);
    void makeTreeIndexExtDFS(int ID, vector<int> &list);
    void ConstructExtensionLabelsNoAllPair();
    void ConstructExtensionLabelsNoAllPairTopDown();
    void ConstructExtensionLabels();
    void ConstructExtensionLabelPartiV(vector<int>& p, bool ifAllPair);
    void ConstructExtensionLabelParti(int pid);
    void ConstructExtensionLabelPartiNoAllPair(int pid);
    void RefreshExtensionLabelsNoAllPair(map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, bool ifParallel);
    void RefreshExtensionLabelsPostMHL(bool ifParallel, bool ifIncrease, double & runT, int threadNum);
    void RefreshExtensionLabels(map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RefreshExtensionLabelPartiV(vector<int>& p, bool ifTopDown);
    void RefreshExtensionLabelParti(int pid);
    void RefreshExtensionLabelPartiTopDown(int pid);
    void RefreshExtensionLabelTopDownV(vector<int>& p);
    void RefreshExtensionLabelTopDown(int ID);
    pair<int,int> ExtensionLabelCompute(int pid, int ID, int ancestor);

    /// VPL
    void VPLGraphPartitionRead(string filename);
    void VPLIndexConstruct();
    void VPLShortcutConstruct();
    void VPLOverlayLabelConstruct();
    void VPLPostLabelConstruct(bool ifParallel, double & t);
    void VPLPartitionIndexConstruct(int pid, bool ifPost);
    void VPLPartitionIndexConstructV(vector<int>& p, bool ifPost);
    void VPLmakeIndexDFSPartiPost(int p, vector<int> &ancestor, vector<int> &interface, map<int, unordered_map<int,int>> &disInfs);
    void VPLCreateIndex_Overlay();
    void VPLmakeIndexDFSOverlay(int p, vector<int>& ancestors);
    void VPLIndexConstructCross(bool ifParallel, double &t);
    void VPLPartitionIndexConstructCross(int pid);
    void VPLPartitionIndexConstructCrossV(vector<int>& p);
    void VPLmakeIndexDFSPartiCross(int p, vector<int>& ancestors);
    int QueryVPL_HASP(int ID1, int ID2);
    int QueryVPL_HASPDebug(int ID1, int ID2);
    int QueryVPL(int ID1, int ID2);
    int QueryVPLDebug(int ID1, int ID2);
    int QueryVPLPartiOverlay(int ID1, int ID2);
    int QueryVPLPartiOverlayDebug(int ID1, int ID2);
    int QueryVPLPartiParti(int ID1, int ID2);
    int QueryVPLSamePartiPost(int ID1, int ID2);
    int QueryVPLSamePartiPostDebug(int ID1, int ID2);
    int QueryVPLHybridCH(int ID1, int ID2);

    void VPLBatchUpdateMixSchedule(vector<pair<pair<int,int>,pair<int,int>>>& wBatchDec, vector<pair<pair<int,int>,pair<int,int>>>& wBatchInc, int batch_i, double &runT1);
    void VPLBatchUpdateMix(vector<pair<pair<int,int>,pair<int,int>>>& wBatchDec, vector<pair<pair<int,int>,pair<int,int>>>& wBatchInc, int batch_i, double &runT1);
    void VPLBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1);
    void VPLShortcutUpdate(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatchDecSchedule, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatchIncSchedule, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchDecSchedule, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchIncSchedule, double& tSc);
    void DecreasePartiBatchVPLShortcut(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC);
    void DecreasePartiBatchVPLShortcutV(vector<int> p, map<int,vector<pair<pair<int,int>,pair<int,int>>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<pair<pair<int,int>,int>>>& updatedSCs);
    void DecreasePartiBatchUpdateCheckVPL(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int,int>,pair<int,int>>>& overlayBatch, bool ifParallel);
    void DecreaseOverlayBatchVPL(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void DecreaseOverlayBatchLabelVPL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL);
    void EachNodeProBDis5VPLOverlay(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);
    void MixRepair_PartiIndexVPLPost(bool ifParallel, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatchDec, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatchInc, double & runT);
    void MixRefreshCrossLabelsVPL(bool ifParallel, double & runT, int threadNum);
    void Repair_PartiIndexVPLPost(bool ifParallel, int updateType, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, double & runT);
    void VPLIndexUpdatePostParti(int pid, bool ifIncrease);
    void VPLIndexUpdatePostPartiV(vector<int>& p, bool ifIncrease);
    void VPLIndexUpdatePostPartiDFS(int pid, int p, vector<int>& ancestor, vector<int> & interface, vector<map<int,int>>& disInfs, set<int>& vertexIDChL, int rootHeight, bool ifIncrease, Timer& tt);
    void RefreshCrossLabelsVPL(bool ifParallel, bool ifIncrease, double & runT, int threadNum);
    void DecreasePartiBatchLabelVPLCrossRoot(int ProBeginVertexID, vector<Node> &Tree, vector<int> &rank, int heightMax, set<int> &vertexIDChL);
    void DecreasePartiBatchLabelVPLCross(int pid, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, set<int> &vertexIDChL);
    void DecreasePartiBatchLabelVPLCrossV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>>& ProBeginVertexSetV, set<int> &vertexIDChL);
    void IncreasePartiBatchLabelVPLCrossRoot(int ProBeginVertexID, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchLabelVPLCross(int pid, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchLabelVPLCrossV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>> &ProBeginVertexSetV, vector<vector<int>> &VidtoTNid);
    void VPLBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1);
    void IncreasePartiBatchUpdateCheckVPL(map<int, vector<pair<pair<int, int>, pair<int, int>>>> &partiBatch, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatch, bool ifParallel);
    void MixVPLBatchShortcutUpdateOverlay(vector<pair<pair<int, int>, pair<int, int>>> &wBatchDec, vector<pair<pair<int, int>, pair<int, int>>> &wBatchInc, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<map<int, vector<pair<int, int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void MixVPLBatchShortcutUpdateParti(map<int, vector<pair<pair<int, int>, pair<int, int>>>> &partiBatchDec, map<int, vector<pair<pair<int, int>, pair<int, int>>>> &partiBatchInc,vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchDec, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchInc, bool ifParallel);
    void MixPartiBatchVPLShortcut(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatchDec, vector<pair<pair<int, int>, pair<int, int>>> &wBatchInc, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSCDec, vector<pair<pair<int, int>, int>> &updatedSCInc);
    void MixPartiBatchVPLShortcutV(vector<int>& p, map<int, pair<vector<pair<pair<int, int>, pair<int, int>>>,vector<pair<pair<int, int>, pair<int, int>>>>>& partiBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<pair<pair<int,int>,int>>>& updatedSCDec, vector<vector<pair<pair<int, int>, int>>> &updatedSCInc);
    void IncreasePartiBatchVPLShortcut(int pid, vector<pair<pair<int, int>, pair<int, int>>> &wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int, int>, int>> &updatedSC);
    void IncreasePartiBatchVPLShortcutV(vector<int> p, map<int, vector<pair<pair<int, int>, pair<int, int>>>> &wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<pair<pair<int, int>, int>>> &updatedSCs);

    void IncreaseOverlayBatchVPL(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<map<int, vector<pair<int, int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void eachNodeProcessIncrease1VPLOverlay(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void IncreaseOverlayBatchLabelVPL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void MixOverlayBatchLabelVPL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void MixVPLPartiLabelUpdateTopDownV(vector<int>& p, vector<map<pair<int,int>,pair<int,int>>>& wBatchDec, vector<map<pair<int, int>, int>> &wBatchInc, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void MixVPLPartiLabelUpdateTopDown(int pid, map<pair<int,int>,pair<int,int>>& wBatchDec, map<pair<int, int>, int> &wBatchInc, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void MixVPLOverlayLabelUpdateTopDown(map<pair<int,int>,pair<int,int>>& wBatchDec, map<pair<int, int>, int> &wBatchInc, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void MixVPLOverlayLabelUpdateBySet(map<pair<int,int>,pair<int,int>>& wBatchDec, map<pair<int, int>, int> &wBatchInc, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void MixVPLPartiLabelUpdateBySetV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void MixVPLPartiLabelUpdateBySet(int pid, map<pair<int,int>,pair<int,int>>& wBatchDec, map<pair<int, int>, int> &wBatchInc, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void TreeNodeMixUpdateVPLOverlay(int children, vector<int>& line, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void TreeNodeMixUpdateVPL(int child, vector<int>& line, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void TreeNodeMixUpdateWhileVPL(int child, vector<int>& line, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void MixPartiBatchCrossLabelVPL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void MixPartiBatchCrossLabelVPLV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void MixPartiBatchLabelVPLCrossRoot(int ProBeginVertexID, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>> &VidtoTNid, set<int> &vertexIDChL);
    void VPLQueryDCHDebug(int ProID, int Cid);
    void VPLQueryCrossDebug(int ProID, int Cid);
    bool judgeQueryAffect(int ID1, int ID2, vector<pair<bool,int>> & partiAffect, int stage_i);

    void VPLBPFlagConstruct(bool ifParallel);
    void BPFlagConstructionParti(int pid);
    void BPFlagConstructionPartiV(vector<int>& p);
    void BPFlagConstructVertex(int ID, unordered_set<int>& targetSet);
    void VPLInnerUnaffectedPartiDetect(vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchDec, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchInc,  bool ifParallel);
    void VPLInnerUnaffectedPartiDetectOverlayDec(vector<bool>& innerAffectTemp);
    void VPLInnerUnaffectedPartiDetectBFS(vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchDec, vector<bool>& innerAffectTemp);
    void VPLOverlayGraphUpdate(vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchDec, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchInc,  bool ifParallel);
    void VPLOverlayUpdateEdgeV(vector<pair<pair<int, int>, pair<int, int>>>& overlayBatch);
    void VPLOverlayUpdateEdge(pair<pair<int, int>, pair<int, int>>& overlayE);
    void VPLOuterUnaffectedPartiDetect();
    void VPLBPFlagUpdate(vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchDec, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchInc, bool ifParallel);
    void BPFlagUpdateIncrease(vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchInc, bool ifParallel);
    void BPFlagUpdateIncParti(int pid);
    void BPFlagUpdateIncPartiV(vector<int>& p);
    void BPFlagUpdateIncVertex(int ID, unordered_set<int>& targetSet, set<pair<int,int>>& PBTemp);
    void BPFlagUpdateDecrease(vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchDec, bool ifParallel);
    void BPFlagUpdateDecreaseNew(vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchDec, bool ifParallel);
    void BPFlagUpdateSPTreeVertex(int ID1, int ID2, vector<int>& pTree);//detect the affected queries by amalgamating SP trees
    void BPFlagUpdateSPTreeVertexTarget(int ID1, int ID2, unordered_set<int> targetSet, vector<int>& pTree);
    void BPFlagUpdateSPTreeEdge(pair<pair<int, int>, pair<int, int>>& overlayE, int& count);
    void BPFlagUpdateSPTreeEdgeSet(pair<pair<int, int>, pair<int, int>>& overlayE, int& count);
    void BPFlagUpdateSPTreeEdgeV(vector<pair<pair<int, int>, pair<int, int>>>& p, int& count);
    void ThreadDistributeBP(vector<pair<pair<int, int>, pair<int, int>>>& overlayBatch, vector<vector<pair<pair<int, int>, pair<int, int>>>>& processID);
    void EllipseSpaceEstimate(int s, int t, int ld, pair<int,int>& resultX, pair<int,int>& resultY);
    void EllipseSpace(int s, int t, int ld, pair<int,int>& resultX, pair<int,int>& resultY);
    pair<double, double> EllipseFunc(double h, double k, double a, double b, double cosA, double sinA, double t);
    void ComputeDynamicMBRForParti(int pid, Rect& result);
    void ComputeOriginalMBRForParti(int pid, Rect& result);
    static bool MySearchCallback(Rect* value);
    void RTreeConstruct();
    void VPLUnaffectedPartiDetectTruth(vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchDecSchedule, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchIncSchedule, double& tDetect, vector<pair<bool,int>>& partiAffectInfo);
    void VPLInnerUnaffectedPartiDetectByOracle(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatchDecSchedule, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatchIncSchedule, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchDecSchedule, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchIncSchedule, int& upNum, double& tDetect, vector<pair<bool,int>>& partiAffectInfo);
    void DistanceOracleTestParti(int pid, vector<bool>& affectedP, unordered_set<int>& affectedPSet, queue<int>& mq);
    void DistanceOracleTestPartiV(vector<int>& p, vector<bool>& affectedP, unordered_set<int>& affectedPSet, queue<int>& mq);
    void EachNodeProBDisUpdateRange(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank, pair<int,int>& p, bool& ProDisChange);
    void EachNodeProBDisUpdate(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank, bool ifParallel);
    void EachNodeProBDisUpdateVert(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank, pair<int,int>& p, bool& ProDisChange);
    void PartitionUpdateClustering(unordered_set<int>& scheduleP, vector<vector<int>>& clusterP);
    void PartitionRegionCluster(int clusterNum, vector<vector<int>>& clusterP, vector<pair<int,int>>& ODpair);
    int PathRetrieval(vector<pair<int,int>>& path);
    void ShortestPathCompare(vector<pair<int,int>>& path1, vector<pair<int,int>>& path2);
    void ReadPathCompare(string file1, string file2);
    void WritePartitionResult();

    /// PostMHL
    void TreeDecompositionPartitioningNaive();
    double TreeDecompositionPartitioning(int pNum, double bRatioUpper, double bRatioLower);
    void TDPContract();
    void TDPCreateTreeAndPartiNaive();
    void TDPCreateTreeAndParti(int pNum, double bRatioUpper, double bRatioLower);
    void MinimizeBoundaryNum(int upperB, int lowerB);
    void GetPartiRootCandidates(int p, int upperB, int lowerB, vector<pair<int,int>>& candidates, map<int,vector<int>>& candidatesMap);
    void GetChildrenNumDFS(int p, vector<int>& ancestors);
    void makePartitionDFS(vector<set<OrderCompMax>>& orderPartiV, int pid, int vid);
    void PostMHLIndexConstructOverlay();
    void PostMHLCreateIndex_Overlay();
    void PostMHLmakeIndexDFSOverlay(int p, vector<int>& ancestors);
    void PostMHLIndexConstructPost(bool ifParallel, double &t);
    void PostMHLPartitionIndexConstructV(vector<int>& p, bool ifPost);
    void PostMHLPartitionIndexConstruct(int pid, bool ifPost);
    void PostMHLmakeIndexDFSPartiPost(int p, vector<int>& ancestor, vector<int> & interface, map<int, unordered_map<int,int>>& disInfs);
    void PostMHLmakeIndexDFSParti(int p, vector<int>& ancestor, vector<int> & interface);
    void PostMHLIndexConstructExtend(bool ifParallel, double &t);
    void PostMHLPartitionIndexConstructExtend(int pid);
    void PostMHLPartitionIndexConstructExtendV(vector<int>& p);
    void PostMHLmakeIndexDFSPartiExtend(int p, vector<int>& ancestors);
    void IndexSizePostMHL();

    void PostMHLIndexStoreCH(string filename);
    void PostMHLIndexCompareCH(string filename);
    void MHLIndexCompareCH(string filename);


    void BoundShortcutsCheck(bool ifParallel, int updateType);
    void BoundShortcutsCheckPartiV(pair<int,int> pRange, vector<int>& pids, vector<vector<pair<pair<int,int>,pair<int,int>>>>& affectedShortcuts, int updateType);
    void BoundShortcutsCheckParti(int pid, vector<pair<pair<int,int>,pair<int,int>>>& affectedShortcut, int updateType);
    void InterfacePropagate(int child, vector<int>& interfaces, vector<Node> &Tree, bool ifIncrease);
    void InterfacePropagateParallel(pair<int,int> pRange, vector<int>& pids, bool ifIncrease);
    void InterfaceEntryDecreaseUpdate(bool ifParallel);
    void AncestorEntryDecreaseUpdateParti(int pid);
    void AncestorEntryDecreaseUpdate(int pid, int child,vector<int>& line, vector<int>& interfaces, set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank, int rootHeight);

    void ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID);
    void ConstructPH2H_PartiV(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_PartiVCH(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_Parti(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
//    void ConstructPH2H_PartiNoLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_PartiLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void H2HCreateTree_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs);//Create tree for partition
    void H2HCreateIndex_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP);//Create labels for partition

    void H2HCreateTree_Overlay();
    void H2HCreateIndex_Overlay();
    int ShortcutDisCheck(int ID1, int ID2);

    void makeTreeIndexDFSP(int p, vector<int>& list,  vector<Node>& TreeP, vector<int>& rankP);

    void makeRMQCoreP(int pid, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs, vector<vector<Node>>& Trees);
    void makeRMQDFSCoreP(int pid, int p, int height, vector<int>& EulerSeqP, vector<vector<int>>& toRMQs, vector<vector<Node>>& Trees);
    void makeRMQCore();
    void makeRMQDFSCore(int p, int height, vector<int>& EulerSeq);
    void makeIndexDFS(int p, vector<int> &list);

    void makeRMQ(vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree);
    void makeRMQDFS(int p, int height, vector<int>& EulerSeq, vector<int>& toRMQ, vector<Node>& Tree);
//    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
//    void insertEMTorder(int u,int v,int w);
    void deleteECore(int u,int v);
    void insertECore(int u,int v,int w);
    void insertECoreMT(int u,int v,int w);
    int matchCore(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank);//vector<pair<int,int>> &vert
    int matchCoreParti(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank);//vector<pair<int,int>> &vert
    void IndexSizePH2H();  //Core-tree index size computation
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void NeighborComorderH2H(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void insertEMTorder(int u,int v,int w);
    void PH2HVertexOrdering(int type);
    void SketchGraphBuild();
    void ReadPartitionGraph();
    void OverlayOrderingBuild();
    void OverlayOrderingBoundaryFirstMDEGreedy();
    void OverlayOrderingBuildBoundaryFirst(int nodenum, vector<vector<int>> Neighbor);
    void PartitionOrderingBuildMDE(bool ifParallel);
    void PartitionOrderingBoundaryFirstMDE(bool ifParallel, vector<bool>& change);
    void OrderingAssemblyMDEBoundaryFirstGreedy(int pNum);
    void OrderingAssemblyMDEBoundaryFirst(int pNum);
    void OrderingAssemblyMDE(int pNum);
    void OrderingAssemblyBoundaryFirst(int pNum);
    void SketchOrder(vector<vector<pair<int,int>>> Neighbor, vector<int> &vNodeOrderSketch);

    void PartitionOrderingV(vector<int>& p);
    void PartitionOrdering(int pid);
    void PartitionOrderingBoundaryFirstV(vector<int>& p, vector<bool>& change);
    void PartitionOrderingBoundaryFirst(int pid, vector<bool>& change);
    void deleteEOrderGenerate(int u,int v);
    void insertEOrderGenerate(int u,int v,int w);



	///Query processing
    int QueryNP(int ID1, int ID2);
    int QueryCHWP(int ID1, int ID2);
    int QueryH2H(int ID1, int ID2);

	//PH2H
    int QueryCore(int ID1, int ID2);
    int QueryH2HPartition(int ID1, int ID2, int PID);
    int QueryH2HPartitionPost(int ID1, int ID2, int PID);
    int QueryPartiCore(int ID1, int ID2);
    int QueryPartiCoreDebug(int ID1, int ID2);
    int QueryPartiCoreExt(int ID1, int ID2);
    int QueryPartiCoreExtLCA(int ID1, int ID2);
    int QuerySameParti(int ID1, int ID2);
    int QuerySamePartiPost(int ID1, int ID2);
    int QuerySamePartiPostOpt(int ID1, int ID2);
    int QueryPartiParti(int ID1, int ID2);
    int QueryPartiPartiExt(int ID1, int ID2);
    int QueryPartiPartiExtLCA(int ID1, int ID2);
    int QueryPartiPartiExtLCADebug(int ID1, int ID2);



    //PCH
    int QueryCoreCH(int ID1, int ID2);
    int QueryCHPartition(int ID1, int ID2, int PID);
    int QueryPartiCoreCH(int ID1, int ID2);
//    int QuerySamePartiCH(int ID1, int ID2);
    int QueryPartiPartiCH(int ID1, int ID2);


    //PostMHL
    int QueryOverlayCH(int ID1, int ID2);
    int QueryOverlayCHDebug(int ID1, int ID2);
    int QueryOverlay(int ID1, int ID2);
    int QueryOverlayDebug(int ID1, int ID2);
    int QueryPostMHLSamePartiPost(int ID1, int ID2);
    int QueryPostMHLSamePartiPostDebug(int ID1, int ID2);
    int QueryPostMHLPartiParti(int ID1, int ID2);
    int QueryPostMHLPartiPartiDebug(int ID1, int ID2);
    int QueryPostMHLPartiPartiExt(int ID1, int ID2);
    int QueryPostMHLPartiPartiExtDebug(int ID1, int ID2);
    int QueryPostMHLPartiOverlay(int ID1, int ID2);
    int QueryPostMHLPartiOverlayDebug(int ID1, int ID2);


    void EffiCheck(string filename,int runtimes);
    void EffiCheckStages(vector<pair<int,int>> & ODpair, int runtimes, int intervalT, unsigned long long & throughputNum, vector<double>& stageUpdateT, vector<double>& stageQueryT);
    vector<double> StageDurationCompute(int intervalT);
    void StageDurationComputeSchedule(int intervalT,vector<vector<double>>& durations);
    void GetBatchThroughput(vector<double> & queryTimes, int intervalT, unsigned long long & throughputNum, vector<double>& stageUpdateT);
    void EffiStageCheck(vector<pair<int,int>> & ODpair, int runtimes, vector<double> & queryTimes);
    void EffiStageCheck(vector<pair<int,int>> & ODpair, int runtimes, vector<vector<double>> & queryTimes, vector<double> & aveQT);
    void EffiStageCheckReal(vector<pair<int,int>> & ODpair, int runtimes, vector<vector<double>> & queryTimes);
    void EffiStageCheck(vector<pair<int,int>> & ODpair, int runtimes, vector<vector<double>> & queryTimes);
    void StagePerformanceShow(int batchNum, vector<double>& stageUpdateT, vector<double>& stageQueryT);
    void AverageStagePerformance(int batchNum, vector<double>& stageUpdateT, vector<double>& stageQueryT);
    double EffiMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType, vector<double> &qTime);
    double EffiPMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType, vector<double> &qTime);
    double EffiPostMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType, vector<double> &qTime);
    double EffiVPLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType, vector<double> &qTime);
    int QueryPostMHL(int ID1, int ID2);
    int QueryPostMHLDebug(int ID1, int ID2);
    int QueryPMHLOpt(int ID1, int ID2);//algoQuery: 0: Dijkstra; 1: PCH-No; 4: PH2H-Post; 5: PH2H-Extend
    int QueryPMHL(int ID1, int ID2);//algoQuery: 0: Dijkstra; 1: PCH-No; 2: PH2H-No; 3: PCH-Post; 4: PH2H-Post; 5: PH2H-Extend
	int Query(int ID1, int ID2);//algoQuery: 0: Dijkstra; 2: PH2H-No; 4: PH2H-Post; 5: PH2H-Extend
    int QueryDebug(int ID1, int ID2);

    int QueryCoreDebug(int ID1, int ID2);
    int LCAQuery(int _p, int _q);
    int LCAQueryPartition(int _p, int _q, int PID);// query within partition
    int LCAQueryPartitionPost(int _p, int _q, int PID);// query within partition
    int LCAQueryOverlay(int _p, int _q);
    int LCAQuery(int _p, int _q, vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree);

    //Correctness Check
    void CorrectnessCheck(int runtimes);
    void CorrectnessCheckCore(int runtimes);
    void DFSTree(vector<int>& tNodes, int id);




    //// For throughput test
    void RealUpdateThroughputTest(string updateFile);
    void RandomUpdateThroughputTest(string updateFile, int batchNum, int batchSize, int batchInterval);
    void RealUpdateThroughputTestQueueModel(string updateFile, int batchNumTest, double T_r, int workerNum, int regionNum);
    void RandomUpdateThroughputTestQueueModel(int batchNum, int batchSize, int batchInterval, double T_r, int workerNum);
    void SPThroughputTest(int updateType, bool ifBatch, int batchNum, int batchSize, int batchInterval, int runtimes);
    void DecBatchThroughput(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);//process batch update
    void IncBatchThroughput(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);//process batch update
    void DecBatchThroughputNP(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);//process batch update
    void EachNodeProBDis5H2H(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis);
    void IncBatchThroughputNP(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);//process batch update
    void eachNodeProcessIncrease1H2H(int children, vector<int>& line, int& changelabel);

    void PMHLBatchUpdateDec(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PMHLBatchUpdateDecOpt(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PostMHLBatchUpdateDec(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PMHLBatchUpdateInc(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PMHLBatchUpdateIncOpt(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PostMHLBatchUpdateInc(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    unsigned long long EffiCheckThroughput(vector<pair<int,int>>& ODpair, Timer& tRecord, int batchInterval, unsigned long long& throughputNum);

    void IndexConstruction();

	/// Index update
    void IndexMaintenance(int updateType, int updateSize, bool ifBatch, int batchSize);
    void DecreaseSingle(int a, int b, int oldW, int newW);//process one update
    void IncreaseSingle(int a, int b, int oldW, int newW);//process one update
    void DecreaseBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//process batch update
    void IncreaseBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//process batch update
    void DecreaseOverlay(int a,int b, int newW, vector<unordered_map<vertex,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void DecreaseParti(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis,
    void EachNodeProBDis5PostMHLOverlay(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis,
    void EachNodeProBDis5Parti(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis, map<int,int>& checkedDis,
    void IncreaseOverlay(int a, int b, int oldW, int newW, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1PostMHL(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1PostMHLOverlay(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void IncreaseParti(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1Parti(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void DecreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);//batch decrease for overlay graph
    void DecreaseH2HBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void DecreaseOverlayBatchPostMHL(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);//batch decrease for overlay graph
    void DecreaseOverlayBatchLabelPostMHL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet,set<int>& vertexIDChL);
    void DecreasePartiBatchLabelPostMHLExtendV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>>& ProBeginVertexSetV, set<int> &vertexIDChL);
    void DecreasePartiBatchLabelPostMHLExtend(int pid, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, set<int> &vertexIDChL);
    void DecreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet,set<int>& vertexIDChL);
    void DecreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay);
    void DecreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC);
    void DecreasePartiBatchUpdateCheckCHV(vector<int> p, map<int,vector<pair<pair<int,int>,pair<int,int>>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<vector<pair<pair<int,int>,int>>>& updatedSC);
    void DecreasePartiBatchUpdateCheckPostMHL(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int,int>,pair<int,int>>>& overlayBatch, bool ifParallel);
    void DecreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void DecreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU, bool ifConstruct);
    void DecreasePartiBatchPostMHLShortcut(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch,  vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC);
    void DecreasePartiBatchPostMHLShortcutV(vector<int> p, map<int,vector<pair<pair<int,int>,pair<int,int>>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<pair<pair<int,int>,int>>>& updatedSCs);
    void DecreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU, bool ifPost);
    void DecreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet,set<int>& vertexIDChL);
    void DecreasePartiBatchLabelV(vector<int> p, vector<vector<Node>> &Trees, vector<vector<int>> &ranks, vector<int> &heightMaxs, vector<vector<int>> &ProBeginVertexSets, vector<set<int>>& vertexIDChLs);
    void IncreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreaseH2HBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreaseOverlayBatchPostMHL(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreaseOverlayBatchLabelPostMHL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay);
    void IncreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt,vector<pair<pair<int,int>,int>>& updatedSC);
    void IncreasePartiBatchUpdateCheckCHV(vector<int> p, map<int,vector<pair<pair<int,int>,pair<int,int>>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& overlayBatch, bool ifOpt, vector<vector<pair<pair<int,int>,int>>>& updatedSCs);
    void IncreasePartiBatchUpdateCheckPostMHL(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatch, bool ifParallel);
    void IncreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void IncreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void IncreasePartiBatchPostMHLShortcut(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch,  vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, vector<int>& PropagateOverlay);
    void IncreasePartiBatchPostMHLShortcutV(vector<int> p, map<int, vector<pair<pair<int, int>, pair<int, int>>>> &wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<pair<pair<int, int>, int>>> &updatedSCs, vector<int>& PropagateOverlay );
    void IncreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU, bool ifPost);
    void IncreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchLabelV(vector<int> p, vector<vector<Node>> &Trees, vector<vector<int>> &ranks, vector<int> &heightMaxs,  vector<vector<int>> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchLabelPostMHLExtendV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>>& ProBeginVertexSetV, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchLabelPostMHLExtend(int pid, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);

	/// Graph Preprocessing
    void ReadGraph(string filename);
    void ReadCoordinate(string filename);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);
    void ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData);
    void ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData);
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
    void ODGeneParti(int num, string filename, int portion);
    void ODGeneSameParti(int num, string filename);
    void ODGeneCrossParti(int num, string filename);
	void UpdateGene(int num, string filename);
    void QueryGenerationParti(bool ifSame);

    int ReadRealUpdates(string filename, vector<vector<pair<pair<int,int>,int>>>& batchUpdates);
    void GetUpdatePortion(string filename);
    void GetUpdatePortion2(int batchInterval);
    void GetQueryDistribution(string filename, int batchInterval);
    vector<int> SampleBatches(string filename, int batchNum);
    void ReadRealQueries(string filename, vector<vector<pair<pair<int,int>,int>>>& ODpair);
    void WriteTreeIndexOverlay(string filename);
    void ReadTreeIndex(string file);
    void WriteTreeIndexParti(string filename);
    void WriteGraph(string graphfile);
    void WriteOrder(string filename);
    void ReadOrder(string filename);
    void CompareOrder(string filename1, string filename2);
    void GraphPartitionRead(string filename);
    void WritePostMHLPartiResult(string filename);
    void ReadPostMHLPartiResult(string filename);

    vector<int> DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);
    vector<int> DFS_CC(vector<vector<pair<int,int>>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);

    void getThroughputEvolveData(string filename);


    float nextTime(double rateParameter) {
        return -log(1.0f - (double) rand() / (RAND_MAX + 1.0)) / rateParameter;
    }
/// poisson process
    std::vector<double> poisson(double lambda, int T) {// T: seconds
        std::vector<double> sequence;
        double t=0;
//    srand(time(NULL));
        int num=0;
        while(t < T) {
            double element = nextTime(lambda);
            if(t+element < T) {
                sequence.push_back(t + element);
                //cout<< t+element<<endl;
            }
            t+=element;
            num++;
        }
        return sequence;
    }
/// function of generating the queries within the duration of T, fulfilling the poisson distribution
    std::vector<query> generate_queryList(double lambda, int T) {
        std::vector<query> list;
        std::vector<double> sequence = poisson(lambda, T);
//        std::cout<<"sequence size of poisson process: "<<sequence.size()<<" ; lambda: "<<lambda<<" ; T: "<<T<<std::endl;
//    std::cout<<"process_time_list size: "<<process_time_list.size()<<std::endl;
//    int m = process_time_list.size();
        int n = sequence.size();
        for(int i=0; i < n; i++) {
            int id = rand()%n;
//        query query(sequence[i], process_time_list[i%m]);
//        list.push_back(query);
            list.emplace_back(sequence[i],0);
//            if(i<2 || i>n-2){
//                std::cout<<i<<": "<<id<<" "<<sequence[i]<<" s"<<std::endl;
//            }
        }
        return list;
    }

    double get_mean(vector<double>& times) {
        double mean = 0.0;
        for (auto val : times) {
            mean += val;
        }
        if(times.size()>0)
            return mean / times.size();
        else
            return 0.0;
    }
    double get_var(vector<double>& times) {
        double mean = get_mean(times);
        double var = 0.0;
        for (auto val : times) {
            var += (val - mean) * (val - mean);
        }
        if(times.size()>0)
            return var / times.size();
        else
            return 0.0;
    }
};

#endif // HEAD_H_
