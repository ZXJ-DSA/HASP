/*
 * Construction.cpp
 *
 *  Created on: 24 August 2023
 *      Author: Xinjie ZHOU
 */
#include "head.h"
#include "PH2H.hpp"

/// Index Construction
void Graph::IndexConstruction(){
    bool flag=false;
    if(threadnum==1){
        cout<<"Single thread computation."<<endl;
        threadnum=40;
        flag=true;
    }
    if(algoChoice==1){
        cout<<"System Index: CH + H2H."<<endl;
        HybridSPIndexConstruct();
    }else if(algoChoice==2){
        cout<<"System Index: PH2H."<<endl;
        HybridPSPIndexConstruct();
//        PH2HIndexConstruct();
    }else if(algoChoice==3){
        cout<<"System Index: PCH + PH2H."<<endl;
        PMHLIndexConstruct();
    }else if(algoChoice==4){
        cout<<"System Index: PCH + PH2H with Optimizations."<<endl;
        PMHLIndexConstructOpt();
    }else if(algoChoice==5){
        cout<<"System Index: PostMHL index."<<endl;
        cout<<"Bandwidth: "<<bandWidth<<endl;
        PostMHLIndexConstruct();
    }else if(algoChoice==6){
        cout<<"System Index: VPL index."<<endl;
        VPLIndexConstruct();
    }else if(algoChoice==0){
        cout<<"A* search."<<endl;
        ReadGraph(sourcePath+"/"+dataset);//
    }
    if(flag==true){
        threadnum=1;
    }
}
//function of hybrid multi-stage SP index construction
void Graph::HybridSPIndexConstruct(){
    string orderfile=sourcePath+"/"+dataset+".order";
//    orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
//    orderfile="/Users/zhouxj/Documents/1-Research/Datasets/beijing/tmp/beijing.PostMHL_100.order";

    double runT1=0, runT2=0, runT3=0;
    Timer tt;
//    ReadGraph(sourcePath+"/"+dataset+".CH6");//
    ReadGraph(sourcePath+"/"+dataset);//

    tt.start();
    MDEContraction(orderfile);
    tt.stop();
    runT1=tt.GetRuntime();
    cout<<"Time for MDE contraction: "<<runT1<<" s."<<endl;

    tt.start();
    makeTree();
    tt.stop();
    runT2=tt.GetRuntime();
    cout<<"Time for Tree construction: "<<runT2<<" s."<<endl;

    tt.start();
    makeIndex();
    tt.stop();
    runT3=tt.GetRuntime();
    cout<<"Time for Index building: "<<runT3<<" s."<<endl;

    cout<<"Overall index construction time: "<<runT1+runT2+runT3<<" s."<<endl;

    IndexsizeH2H();
}

//function for MDE contraction
void Graph::MDEContraction(string orderfile){
    cout<<"MDE contraction..."<<endl;
    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }
    ifstream IF(orderfile);
//    if(true){
    if(!IF.is_open()){/// if no order file, use MDE to generate order
        cout<<"Cannot open vertex ordering file "<<orderfile<<endl;
        int Twidth=0;//tree width
        //initialize SCconNodesMT
        SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());

        //initialize E
        map<int,pair<int,int>> m;
        E.assign(node_num,m);
        for(int i=0;i<Neighbor.size();i++){
            for(int j=0;j<Neighbor[i].size();j++)
                E[i].insert(make_pair(Neighbor[i][j].first,make_pair(0,1)));
        }

        _DD_.assign(node_num,0);
        DD.assign(node_num,0);

        set<DegComp> Deg;
        int degree;
        for(int i=0;i<node_num;i++){
            degree=Neighbor[i].size();
            if(degree!=0){
                _DD_[i]=degree;
                DD[i]=degree;
                Deg.insert(DegComp(i));
            }
        }

        vector<bool> exist; exist.assign(node_num,true);
        vector<bool> change; change.assign(node_num,false);

        vector<pair<int,pair<int,int>>> vect;
        NeighborCon.assign(node_num,vect); //NeighborCon.clear();
        //SCconNodes.clear();

        //cout<<"Begin to contract"<<endl;
        int count=0;

        while(!Deg.empty()){
            if(count%100000==0)
                cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
            count+=1;
            int x=(*Deg.begin()).x;

            while(true){
                if(change[x]){
                    Deg.erase(DegComp(x));
                    _DD_[x]=DD[x];
                    Deg.insert(DegComp(x));
                    change[x]=false;
                    x=(*Deg.begin()).x;
                }else
                    break;
            }

            vNodeOrder.push_back(x);
            Deg.erase(Deg.begin());
            exist[x]=false;

            vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

            for(auto it=E[x].begin();it!=E[x].end();it++){
                if(exist[(*it).first]){
                    Neigh.push_back(*it);
                }
            }

            if(Neigh.size()>Twidth)
                Twidth=Neigh.size();

            NeighborCon[x].assign(Neigh.begin(),Neigh.end());

            //multi threads for n^2 combination
            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteEOrderGenerate(x,y);
                change[y]=true;
            }

            int stepf=Neigh.size()/threadnum;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*stepf;
                if(i==threadnum-1)
                    p.second=Neigh.size();
                else
                    p.second=(i+1)*stepf;
                threadf.add_thread(new boost::thread(&Graph::NeighborComOrderGenerate, this, boost::ref(Neigh), p));
            }
            threadf.join_all();
        }

        NodeOrder.assign(node_num,-1);
        for(int k=0;k<vNodeOrder.size();k++){
            NodeOrder[vNodeOrder[k]]=k;
        }
        ofstream ofile(orderfile);
        ofile << node_num << endl;
        for(int i = 0; i < NodeOrder.size(); i++)
            ofile << i << " " << NodeOrder[i] << endl;
        ofile.close();
        cout<<"Finish Contract"<<" , treewidth "<<Twidth<<endl;
        exit(0);
    }
    else{///if there is an order file
        cout<<"Reading vertex ordering... "<< orderfile<<endl;
        NodeOrder.assign(node_num, -1);
        vNodeOrder.assign(node_num, -1);
        int num, nodeID, nodeorder;
        string line;
        getline(IF,line);
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        num=stoi(vs[0]);
        for(int i=0;i<num;i++){
            vs.clear();
            getline(IF,line);
            boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
            if(vs.size()!=2){
                cout<<"Wrong syntax for ordering. "<<line<<endl; exit(1);
            }
            nodeID=stoi(vs[0]);nodeorder=stoi(vs[1]);

            NodeOrder[nodeID]=nodeorder;
            if(nodeorder!=-1){
                vNodeOrder[nodeorder]=nodeID;
            }else{
                cout<<"Wrong order! "<<nodeID<<" "<<nodeorder<<endl; exit(1);
            }
        }
        IF.close();
        unordered_set<int> vertices; vertices.clear();
        for(int i=0;i<node_num;++i){
            if(vertices.find(vNodeOrder[i])==vertices.end()){//if not found
                vertices.insert(vNodeOrder[i]);
            }
        }
        if(vertices.size()!=node_num){
            cout<<"Order wrong! "<<vertices.size()<<" "<<node_num<<endl; exit(1);
        }
        vertices.clear();

        for(int i=0;i<2;++i){
            int id=vNodeOrder[node_num-1-i];
            cout<<"Order "<<node_num-1-i<<": "<<id<<" "<<Neighbor[id].size()<<endl;
        }
        for(int i=0;i<2;++i){
            int id=vNodeOrder[i];
            cout<<"Order "<<i<<": "<<id<<" "<<Neighbor[id].size()<<endl;
        }

        vector<pair<int,pair<int,int>>> vect;
        NeighborCon.assign(node_num,vect);

        SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());//record the supportive vertices of a shortcut, only record edge once by leveraging the ID positions of endpoints

        //initialize E
        map<int,pair<int,int>> m;
        E.assign(node_num,m);
        for(int i=0;i<Neighbor.size();i++){
            for(int j=0;j<Neighbor[i].size();j++)
                E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
        }

        vector<bool> exist; exist.assign(node_num,true);
        //vector<bool> change; change.assign(nodenum,false);

        //cout<<"Begin to contract"<<endl;
        for(int nodeorder=0;nodeorder<node_num;nodeorder++){//start from the most important vertex
            int x=vNodeOrder[nodeorder];
            if(x!=-1){//to identify and exclude the isolated vertices
                exist[x]=false;

                vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

                for(auto it=E[x].begin();it!=E[x].end();it++){
                    if(exist[(*it).first]){
                        Neigh.push_back(*it);
                    }
                }
                NeighborCon[x].assign(Neigh.begin(),Neigh.end());

                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteEorder(x,y);
                    //change[y]=true;
                }

                int ID1,ID2;
                if(Neigh.size()<=100){
//                if(true){
                    //single thread
                    for(int i=0;i<Neigh.size();i++){
                        ID1=Neigh[i].first;
                        for(int j=i+1;j<Neigh.size();j++){
                            ID2=Neigh[j].first;
                            int temp=Neigh[i].second.first+Neigh[j].second.first;
                            insertEorder(ID1,ID2,temp);
                            if(ID1<ID2){
                                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                                    SCconNodesMT[ID1].insert({ID2,vector<pair<int,int>>()});
                                }
                                SCconNodesMT[ID1][ID2].emplace_back(x,temp);//only record onece
                            }
                            else if(ID2<ID1){
                                if(SCconNodesMT[ID2].find(ID1)==SCconNodesMT[ID2].end()){//if not found
                                    SCconNodesMT[ID2].insert({ID1,vector<pair<int,int>>()});
                                }
                                SCconNodesMT[ID2][ID1].emplace_back(x,temp);
                            }

                        }
                    }
                }else{
                    if(Neigh.size()>threadnum){
                        int step=Neigh.size()/threadnum;
                        boost::thread_group thread;
                        for(int i=0;i<threadnum;i++){
                            pair<int,int> p;
                            p.first=i*step;
                            if(i==threadnum-1)
                                p.second=Neigh.size();
                            else
                                p.second=(i+1)*step;
                            thread.add_thread(new boost::thread(&Graph::NeighborComorderH2H, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }else{
                        boost::thread_group thread;
                        for(int i=0;i<Neigh.size();i++){
                            pair<int,int> p;
                            p.first=i; p.second=(i+1);
                            thread.add_thread(new boost::thread(&Graph::NeighborComorderH2H, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }

                }

            }
            else{
                cout<<"Wrong order! "<<x<<" "<<nodeorder<<endl; exit(1);
            }
        }
    }
    NodeOrder_ = NodeOrder;
}

void Graph::insertEMTorder(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }
}

void Graph::NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
//            vSm[ID1]->wait();
//            vSm[ID2]->wait();
            insertEMTOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
//            vSm[ID1]->notify();
//            vSm[ID2]->notify();
        }
    }
//    sm->notify();
}

void Graph::insertEMTOrderGenerate(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//        DD2[u]++;
    }
}

void Graph::deleteEorder(int u,int v){
//    if(E[u].find(v)!=E[u].end()){
//        E[u].erase(E[u].find(v));
//        //DD[u]--;
//    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        //DD[v]--;
    }
}

void Graph::insertEorder(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        //DD[u]++;
        //DD2[u]++;
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        //DD[v]++;
        //DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}

//function for H2H tree construction
void Graph::makeTree(){
    cout<<"Building H2H tree..."<<endl;
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(node_num,vecemp);//record the tree node id whose unique vertex involves this vertex as neighbor

    rank.assign(node_num,0);
    //Tree.clear();
    int len=vNodeOrder.size()-1;
    heightMax=0;

    Node rootn;
    int x=vNodeOrder[len];//from the highest vertex
    //cout<<"len "<<len<<" , ID "<<x<<endl;
    while(x==-1){//to skip those vertices whose ID is -1
        len--;
        x=vNodeOrder[len];
        //cout<<"len "<<len<<" , ID "<<x<<endl;
    }
    rootn.vert=NeighborCon[x];
    rootn.uniqueVertex=x;
    rootn.pa=-1;
    rootn.height=1;
    rank[x]=0;
    Tree.push_back(rootn);
    len--;

    int nn;
    for(;len>=0;len--){
        int x=vNodeOrder[len];
        Node nod;
        nod.vert=NeighborCon[x];
        nod.uniqueVertex=x;
        int pa=match(x,NeighborCon[x]);
        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        nod.hdepth=Tree[pa].height+1;

        for(int i=0;i<NeighborCon[x].size();i++){
            nn=NeighborCon[x][i].first;
            VidtoTNid[nn].push_back(Tree.size());
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1)
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
        }
        if(nod.height>heightMax) {
            heightMax=nod.height;
        }
        rank[x]=Tree.size();

        Tree.push_back(nod);
        //cout<<"len "<<len<<" , ID "<<x<<endl;
    }
}
int Graph::match(int x,vector<pair<int,pair<int,int>>> &vert){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(rank[vert[i].first]>rank[nearest])
            nearest=vert[i].first;
    }
    int p=rank[nearest];
    return p;
}
//function of H2H index construction
void Graph::makeIndex(){
    cout<<"Building H2H index..."<<endl;
    makeRMQ(toRMQ, RMQIndex, Tree);//build LCA index

    //initialize
    vector<int> list; //list.clear();
    list.push_back(Tree[0].uniqueVertex);
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);

    for(int i=0;i<Tree[0].ch.size();i++){
        makeIndexDFS(Tree[0].ch[i],list);
    }

}

//function for computing the index size
void Graph::IndexsizeH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;
    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*2*sizeof(int);//dis
        m3+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].cnt.size()*sizeof(int);//cnt
        m5+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count

    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4+m5;
    cout<<"CH index size: "<<(double)(m4+m5)/1024/1024<<" MB"<<endl;
    cout<<"CH Update information size: "<<(double)m4/1024/1024<<" MB"<<endl;
    cout<<"Distance label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"H2H label size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"H2H Update information size: "<<(double)(m2+m4)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::PMHLIndexConstruct() {
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
//    orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDEGreedy";
    ReadOrder(orderfile);

    string partitionfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
    tt.start();
    /// Partition index and Overlay graph construction
//    Construct_PartiIndex(false);
    Construct_PartiIndex(true, true);
//    PMHLConstructPartiIndexCH(true);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition index construction time: "<<runT1<<" s"<<endl;

    /// Overlay graph construction
    tt.start();
//    Construct_OverlayGraph(true);
    Construct_OverlayGraphNoAllPair(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
    Construct_OverlayIndex(true);
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


    /// Partition index repair
    if(algoUpdate>=PH2H_Post){
        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true, true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
//        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)

    if(algoUpdate==PH2H_Cross){
        tt.start();
//        ConstructExtensionLabelsNoAllPair();//correct version
        ConstructExtensionLabelsNoAllPairTopDown();
        tt.stop();
        runT5=tt.GetRuntime();
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
    }

    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}
//Index construction for VPL
void Graph::VPLIndexConstruct() {
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
//    orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDEGreedy";
    ifstream inFile(orderfile, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << orderfile << endl;
        inFile.close();
        PH2HVertexOrdering(3);//Boundary-first heuristic ordering
        exit(0);
    }else{
        inFile.close();
        ReadOrder(orderfile);
    }
    if(dataset!="Test"){
        ReadCoordinate(sourcePath+"/"+dataset+".co");
    }

    string partitionfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
//    partitionfile=sourcePath+"/tmp";
    VPLGraphPartitionRead(partitionfile);//read partitions
//    GetUpdatePortion(sourcePath+"/"+dataset+"_20160105_300.batchUpdates");
//    GetUpdatePortion2(300);
//    GetQueryDistribution(sourcePath+"/"+dataset+".realQueries", 300);

    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }

    Timer tt;
    /// Shortcut tree construction
    tt.start();
    VPLShortcutConstruct();
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Shortcut tree construction time: "<<runT1<<" s"<<endl;

    tt.start();
//    VPLBPFlagConstruct(true);
//    VPLBPFlagConstruct(false);
    innerAffectedParti.assign(partiNum,true);
    AffectedParti.assign(partiNum,true);
    tt.stop();
    runT2 = tt.GetRuntime();
//    cout<<"BP-Flag construction time: "<<runT2<<" s"<<endl;

    /// Overlay index construction
    tt.start();
    VPLOverlayLabelConstruct();
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


    /// Partition index construction, single thread
    if(algoUpdate>=PH2H_Post){
        double tPost=0;
        tt.start();
        VPLPostLabelConstruct(true,tPost);
//        VPLPostLabelConstruct(false,tPost);
        tt.stop();
        runT4 = tt.GetRuntime();
        algoQuery=PH2H_Post;
        cout<<"Post-boundary label construction time: "<<tt.GetRuntime()<<" s."<<endl;
    }
    if(algoUpdate==PH2H_Cross){
        double tCross=0;
        tt.start();
        VPLIndexConstructCross(true,tCross);
//        VPLIndexConstructCross(false,tCross);
        tt.stop();
        runT5=tt.GetRuntime();
        algoQuery=PH2H_Cross;
        cout<<"Cross-boundary label construction time: "<<runT5<<" s"<<endl;
    }
    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;
    /// print the tree
    /*queue<int> q;
//    stack<int> q;
    q.push(Tree[0].uniqueVertex);
    int ID;
    int h=0;
    while(!q.empty()){
        ID=q.front();
//        ID=q.top();
        q.pop();
        if(Tree[rank[ID]].height==h){
            cout<<ID+1<<"(pos:";
            for(int i=0;i<Tree[rank[ID]].pos.size();++i){
                cout<<Tree[rank[ID]].pos[i]<<",";
            }
            cout<<"dis:";
            for(int i=0;i<Tree[rank[ID]].dis.size();++i){
                cout<<Tree[rank[ID]].dis[i]<<",";
            }
            cout<<"anc:";
            for(int i=0;i<Tree[rank[ID]].vAncestor.size();++i){
                cout<<Tree[rank[ID]].vAncestor[i]+1<<",";
            }
            cout<<"N:";
            for(int i=0;i<Tree[rank[ID]].vert.size();++i){
                cout<<Tree[rank[ID]].vert[i].first+1<<",";
            }
            cout<<"sc:";
            for(int i=0;i<Tree[rank[ID]].vert.size();++i){
                cout<<Tree[rank[ID]].vert[i].second.first<<",";
            }
            cout<<") ";
        }else{
            h=Tree[rank[ID]].height;
            cout<<endl;
            cout<<h<<": "<<ID+1<<"(pos:";
            for(int i=0;i<Tree[rank[ID]].pos.size();++i){
                cout<<Tree[rank[ID]].pos[i]<<",";
            }
            cout<<"dis:";
            for(int i=0;i<Tree[rank[ID]].dis.size();++i){
                cout<<Tree[rank[ID]].dis[i]<<",";
            }
            cout<<"anc:";
            for(int i=0;i<Tree[rank[ID]].vAncestor.size();++i){
                cout<<Tree[rank[ID]].vAncestor[i]+1<<",";
            }
            cout<<"N:";
            for(int i=0;i<Tree[rank[ID]].vert.size();++i){
                cout<<Tree[rank[ID]].vert[i].first+1<<",";
            }
            cout<<"sc:";
            for(int i=0;i<Tree[rank[ID]].vert.size();++i){
                cout<<Tree[rank[ID]].vert[i].second.first<<",";
            }
            cout<<") ";
        }
        for(int i=0;i<Tree[rank[ID]].ch.size();++i){
            q.push(Tree[Tree[rank[ID]].ch[i]].uniqueVertex);
        }
    }
    cout<<endl;*/
    IndexSizePostMHL();//index (labeling+pruning point) size computation
    /// R-Tree computation
    RTreeConstruct();
}

void Graph::PMHLIndexConstructOpt() {
    double runT0, runT1, runT2, runT3, runT4, runT5;
    runT0=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
    orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    ReadOrder(orderfile);

    string partitionfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
//    tt.start();
//    PreConstructAllPairs(true);
//    tt.stop();
//    runT0 = tt.GetRuntime();
//    cout<<"All-pair boundary edges insertion time: "<<runT0<<" s"<<endl;

    tt.start();
    /// Partition index and Overlay graph construction
    Construct_PartiIndex(true, false);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition shortcuts construction time: "<<runT1<<" s"<<endl;

    /// Overlay graph construction
    tt.start();
    Construct_OverlayGraphNoAllPair(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
    Construct_OverlayIndex(true);
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay shortcuts construction time: "<<runT3<<" s."<< endl;
    algoQuery=PCH_No;

    /// Partition index repair
    if(algoUpdate>=PH2H_Post){
        /// partition shortcuts update and label construction
        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true, true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;
        algoQuery=PH2H_Post;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)

    if(algoUpdate==PH2H_Cross){
        tt.start();
//        ConstructExtensionLabelsNoAllPair();//correct version
        ConstructExtensionLabelsNoAllPairTopDown();
        tt.stop();
        runT5=tt.GetRuntime();
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
        algoQuery=PH2H_Cross;
    }

    cout<<"Overall Construction Time: "<<runT0+runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}

void Graph::PostMHLIndexConstruct() {
    double runT0, runT1, runT11, runT2, runT3, runT4, runT5;
    runT0=0, runT1=0, runT11=0, runT2=0, runT3=0, runT4=0, runT5=0;
    Timer tt;
    ReadGraph(sourcePath+"/"+dataset);//
//    ReadGraph(sourcePath+"/"+dataset+".CH19");//
//    ReadGraph(sourcePath+"/"+dataset+".CH2");//
//    ReadGraph(sourcePath+"/"+dataset+".time");//

    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }
//    tt.start();
//    PreConstructAllPairs(true);
//    tt.stop();
//    runT0 = tt.GetRuntime();
//    cout<<"All-pair boundary edges insertion time: "<<runT0<<" s"<<endl;

    tt.start();
    /// Partition and Overlay graph construction
//    TreeDecompositionPartitioningNaive();
    runT11=TreeDecompositionPartitioning(partiNum,bRatioUpper,bRatioLower);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Tree Decomposition time: "<<runT11<<" s"<<endl;
    cout<<"Tree-Decomposition based Partitioning time: "<<runT1<<" s"<<endl;
    // write the partition results to disk
    string partiF = sourcePath+"/tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".partiInfo";
    ifstream IF(partiF);
    if(!IF){//if the label file does not exist, construct it
        WritePostMHLPartiResult(sourcePath+"/tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".partiInfo");
    }
    IF.close();

    algoQuery=PCH_No;


    /// Overlay index construction
    tt.start();
    PostMHLIndexConstructOverlay();
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;
//    algoQuery=PCH_No;

    /// Partition index construction, single thread
    if(algoUpdate>=PH2H_Post){
        double tPost=0;
        tt.start();
        PostMHLIndexConstructPost(true,tPost);
//        PostMHLIndexConstructPost(false);
        tt.stop();
        runT4 = tt.GetRuntime();
        algoQuery=PH2H_Post;
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;
    }
    if(algoUpdate==PH2H_Cross){
        double tCross=0;
        tt.start();
        PostMHLIndexConstructExtend(true,tCross);
//        PostMHLIndexConstructExtend(false);
        tt.stop();
        runT5=tt.GetRuntime();
        algoQuery=PH2H_Cross;
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
    }


    cout<<"Construction Time with partitioning: "<<runT0+runT1+runT2+runT3+runT4+runT5<<" s."<<endl;
    cout<<"Overall Construction Time: "<<runT0+runT11+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePostMHL();
}



// partitioned SP
void Graph::HybridPSPIndexConstruct(){
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
    orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    ReadOrder(orderfile);

    string partitionfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
    tt.start();
    /// Partition index and Overlay graph construction
    Construct_PartiIndex(true, true);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition index construction time: "<<runT1<<" s"<<endl;


    /// Overlay graph construction
    tt.start();
    Construct_OverlayGraph(true);
//    Construct_OverlayGraphNoAllPair(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
    Construct_OverlayIndex(true);
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


    /// Partition index repair
    if(algoUpdate>=PH2H_Post){
        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true, true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)

    if(algoUpdate==PH2H_Cross){
        tt.start();
        ConstructExtensionLabels();
        tt.stop();
        runT5=tt.GetRuntime();
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
    }

    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}

void Graph::PH2HIndexConstruct(){
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
    orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    ReadOrder(orderfile);

    string partitionfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
    tt.start();
    /// Partition index and Overlay graph construction
    Construct_PartiIndex(true,true);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition index construction time: "<<runT1<<" s"<<endl;


    /// Overlay graph construction
    tt.start();
    Construct_OverlayGraph(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
    Construct_OverlayIndex(true);
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


    /// Partition index repair
    if(algoUpdate>=PH2H_Post){
        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true, true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)

    if(algoUpdate==PH2H_Cross){
        tt.start();
        ConstructExtensionLabels();
        tt.stop();
        runT5=tt.GetRuntime();
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
    }

    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}


//function for computing the index size
void Graph::IndexSizePH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;

    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*sizeof(int);//dis
        m1+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].dis.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count
    }
    //overlay
    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m2+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    if(algoUpdate>=PCH_No){
        for(int pid=0;pid<Trees.size();++pid){
            for(int i=0;i<Trees[pid].size();i++){
//                    m3+=Trees[pid][i].dis.size()*sizeof(int);//dis
//                    m3+=Trees[pid][i].pos.size()*sizeof(int);//pos
//                    m4+=Trees[pid][i].dis.size()*sizeof(int);//cnt
                m4+=Trees[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
            }
        }
    }

    if(algoUpdate>=PH2H_No){
        //partitions
        for(int pid=0;pid<Trees.size();++pid){
            for(int i=0;i<Trees[pid].size();i++){
                m3+=Trees[pid][i].dis.size()*sizeof(int);//dis
                m3+=Trees[pid][i].pos.size()*sizeof(int);//pos
                m4+=Trees[pid][i].dis.size()*sizeof(int);//cnt
//                m4+=Trees[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
            }
        }
        //partitions
        for(int i=0;i< SCconNodesMTP.size();i++){
            for(auto it=SCconNodesMTP[i].begin(); it!=SCconNodesMTP[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }
    }

    if(algoUpdate>=PH2H_Post){
        for(int pid=0;pid<TreesPost.size();++pid){
            for(int i=0;i<TreesPost[pid].size();i++){
                m3+=TreesPost[pid][i].dis.size()*sizeof(int);//dis
                m3+=TreesPost[pid][i].pos.size()*sizeof(int);//pos
                m4+=TreesPost[pid][i].dis.size()*sizeof(int);//cnt
                m4+=TreesPost[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
            }
        }
        for(int i=0;i< SCconNodesMTPost.size();i++){
            for(auto it=SCconNodesMTPost[i].begin(); it!=SCconNodesMTPost[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }
    }
    if(algoUpdate==PH2H_Cross){
        for(int i=0;i<TreeExt.size();i++){
            m5+=TreeExt[i].dis.size()*sizeof(int);//dis
            m5+=TreeExt[i].pos.size()*sizeof(int);//pos
            m5+=TreeExt[i].cnt.size()*sizeof(int);//cnt
        }

    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4+m5;
    cout<<"Distance labeling size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"Overlay graph index size: "<<(double)(m1+m2)/1024/1024<<" MB"<<endl;
    cout<<"Partition graphs index size "<<(double)(m3+m4)/1024/1024<<" MB"<<endl;
    cout<<"Extended label size "<<(double)(m5)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

//function for computing the index size
void Graph::IndexSizePostMHL(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;

    //index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*2*sizeof(int);//dis
        m1+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].cnt.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count
        m3+=Tree[i].disInf.size()*2*sizeof(int);
        m3+=Tree[i].disPost.size()*2*sizeof(int);
    }
    //
    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m2+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }


    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4+m5;
    cout<<"Distance labeling size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::PH2HVertexOrdering(int type){
//    ReadGraph(sourcePath+"/"+dataset);
    string filename=sourcePath+"/"+dataset;
    ifstream inGraph(filename);
    if(!inGraph){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    cout<<"Reading graph... "<<filename<<endl;
    Timer tt;
    tt.start();
    string line;
    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    node_num=stoi(vs[0]); edge_num=stoull(vs[1]);
    inGraph.close();
    int pNum=partiNum;
    switch (type) {
        case 0:{//MDE partition + distant MDE overlay
            cout<<"MDE ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuild();
            PartitionOrderingBuildMDE(true);
//            PartitionOrderingBuildMDE(false);
            OrderingAssemblyMDE(pNum);
            break;
        }
        case 1:{
            cout<<"Boundary-first ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuildBoundaryFirst(node_num, NeighborSketch);
            PartitionOrderingBuildMDE(true);
            OrderingAssemblyBoundaryFirst(pNum);
            break;
        }
        case 2:{
            cout<<"Boundary-first MDE ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuild();
            PartitionOrderingBuildMDE(true);
            OrderingAssemblyMDEBoundaryFirst(pNum);
            break;
        }
        case 3:{
            cout<<"Greedy Boundary-first MDE ordering."<<endl;
            ReadPartitionGraph();
            vector<bool> change; change.assign(node_num,false);
            PartitionOrderingBoundaryFirstMDE(true, change);
//            PartitionOrderingBoundaryFirstMDE(false, change);
            OverlayOrderingBoundaryFirstMDEGreedy();
            OrderingAssemblyMDEBoundaryFirstGreedy(pNum);
            break;
        }
        default:{
            cout<<"Wrong ordering type! "<<type<<endl; exit(1);
        }
    }
    tt.stop();
    cout<<"Overall ordering time: "<<tt.GetRuntime()<<" s."<<endl;
    exit(0);
}
void Graph::OrderingAssemblyMDEBoundaryFirstGreedy(int pNum){
    string filename=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDEGreedy";

    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int ID;
    NodeOrder.assign(node_num,-1);
    vNodeOrder.clear();
    int order_i=0;
    set<int> vertices;
    /// For non-boundary vertex
    for(int pid=0;pid<partiNum;++pid){
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[pid].begin();it!=vNodeOrderParti[pid].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            if(!PartiTag[ID].second){// if not boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }else{
                cout<<"Wrong. "<<ID<<endl;
            }
        }
    }
    /// For boundary vertex
    for(auto it=vNodeOrderOverlay.begin();it!=vNodeOrderOverlay.end();++it){
        ID=*it;
        if(PartiTag[ID].second){// if boundary
            vNodeOrder.emplace_back(ID);
            vertices.insert(ID);
            ++order_i;
        }else{
            cout<<"Wrong. "<<ID<<endl;
        }
    }


//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num || vertices.size()!=node_num || vNodeOrder.size()!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<vertices.size()<<" "<<vNodeOrder.size()<<" "<<node_num<<endl; exit(1);
    }



    for(int i=0;i<node_num;++i){
        NodeOrder[vNodeOrder[i]]=i;
    }

    for(int i=0;i<nodeNumOverlay;++i){
        ID=vNodeOrder[node_num-i-1];
        if(!PartiTag[ID].second){
            cout<<nodeNumOverlay<<", "<<node_num-i<<": "<<ID<<"("<<NodeOrder[ID]<<") is not boundary vertex!"<<endl; exit(1);
        }
    }

    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}

void Graph::OrderingAssemblyMDEBoundaryFirst(int pNum){
    string filename=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE2";

    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int PID,ID;
    NodeOrder.assign(node_num,-1);
    vNodeOrder.clear();
    int order_i=0;
    set<int> vertices;
    /// For non-boundary vertex
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            if(!PartiTag[ID].second){// if not boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }
        }
    }
    /// For boundary vertex
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            if(PartiTag[ID].second){// if boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }
        }
    }

//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num || vertices.size()!=node_num || vNodeOrder.size()!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<vertices.size()<<" "<<vNodeOrder.size()<<" "<<node_num<<endl; exit(1);
    }



    for(int i=0;i<node_num;++i){
        NodeOrder[vNodeOrder[i]]=i;
    }

    for(int i=0;i<nodeNumOverlay;++i){
        ID=vNodeOrder[node_num-i-1];
        if(!PartiTag[ID].second){
            cout<<nodeNumOverlay<<", "<<node_num-i<<": "<<ID<<"("<<NodeOrder[ID]<<") is not boundary vertex!"<<endl; exit(1);
        }
    }

    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}
//function of MDE ordering assemblying
void Graph::OrderingAssemblyMDE(int pNum){
    string filename=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE";

    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int PID,ID;
    NodeOrder.assign(node_num,-1);
    int order_i=0;
    set<int> vertices;
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            vertices.insert(ID);
            NodeOrder[ID]=order_i;
            ++order_i;
        }
    }
//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<node_num<<endl; exit(1);
    }
    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}
//function of boundary-first assemblying
void Graph::OrderingAssemblyBoundaryFirst(int pNum){
    string orderfile=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(pNum)+"/vertex_order2";
    set<int> vcheck;//use to check the redundant ordered vertex
    vcheck.clear();
    vNodeOrder.clear();
    int pid,ID;
    //for vertex within partition
    for(int k=0;k<vNodeOrderOverlay.size();k++){
        pid=vNodeOrderOverlay[k];

        for(auto it=vNodeOrderParti[k].begin();it!=vNodeOrderParti[k].end();++it){
//            ID=it->second;
            ID=*it;
            if(!PartiTag[ID].second){//if not boundary vertex
                vNodeOrder.push_back(ID);
                if(vcheck.find(ID)!=vcheck.end())
                    cout<<"wrong: redundant vertex ordered"<<endl;
                vcheck.insert(ID);
            }
        }
    }

    //for boundary vertex
    for(int k=0;k<vNodeOrderOverlay.size();k++){
        pid=vNodeOrderOverlay[k];
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID=BoundVertex[pid][i];
            vNodeOrder.push_back(ID);
            if(vcheck.find(ID)!=vcheck.end())
                cout<<"wrong: redundant vertex ordered"<<endl;
            vcheck.insert(ID);
        }
    }

    //cout<<"total number of ordered vertex "<<vNodeOrder.size()<<endl;
    if(vNodeOrder.size()!=node_num)
        cout<<"Something wrong happened: some vertices do not have the vertex order!"<<endl;

    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }

    ofstream OF(orderfile);
    if(!OF){
        cout<<"Cannot open file "<<orderfile<<endl;
        exit(1);
    }
    OF<<NodeOrder.size()<<endl;
    for(int i=0;i<NodeOrder.size();i++){
        OF<<i<<" "<<NodeOrder[i]<<endl;
    }
    OF.close();
    cout<<"Finished."<<endl;
}

void Graph::SketchOrder(vector<vector<pair<int,int>>> Neighbor, vector<int> &vNodeOrderSketch){
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(1,1)));
    }

    _DD_.assign(node_num,0);
    DD.assign(node_num,0);

    set<DegComp> Deg;
    int degree;
    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();
        if(degree!=0){
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
        }
    }

    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect); //NeighborCon.clear();
    //SCconNodes.clear();

    //cout<<"Begin to contract"<<endl;
    int count=0;

    while(!Deg.empty()){
        //if(count%10==0) cout<<"count "<<count<<endl;
        count+=1;
        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }

        vNodeOrderSketch.push_back(x);
        Deg.erase(Deg.begin());
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();i++){
            for(int j=i+1;j<Neigh.size();j++){
                insertEOrderGenerate(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                change[Neigh[i].first]=true;
                change[Neigh[j].first]=true;
            }
        }
    }

    /*NodeOrderSketch.assign(nodenum,-1);
    for(int k=0;k<vNodeOrderSketch.size();k++){
        NodeOrderSketch[vNodeOrderSketch[k]]=k;
        cout<<"order "<<k<<", vertex "<<vNodeOrderSketch[k]<<", degree "<<NeighborSketch[vNodeOrderSketch[k]].size()<<endl;
    }*/
    //cout<<"Finish Contract"<<endl;
}

void Graph::SketchGraphBuild(){
    NeighborsParti.assign(node_num, vector<pair<vertex,int>>());
    NeighborsOverlay.assign(node_num,unordered_map<vertex,int>());
//    NeighborsOverlay.assign(node_num,vector<pair<vertex,int>>());
    PartiTag.assign(node_num, make_pair(-1,false));

    bool flag_minus = false;

    string filename=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    ifstream IF1(filename+"/subgraph_vertex");
    if(!IF1){
        cout<<"Cannot open file "<<"subgraph_vertex"<<endl;
        exit(1);
    }

    int pnum2;
    IF1>>pnum2;
    if(algoParti == "NC"){
//        flag_minus = true;
        partiNum = pnum2;
    }else if(algoParti == "SC" || algoParti == "MT"){
//        flag_minus = true;
//        pnum2 = pnum;
    }
    cout<<"Partition number: "<<pnum2<<endl;

    PartiVertex.assign(partiNum,vector<vertex>());
    for(int k=0;k<pnum2;k++){
        int vernum,ID;
        IF1>>vernum;
        for(int i=0;i<vernum;i++){
            IF1>>ID;
//            if(flag_minus){
//                ID = ID-1;
//            }

            if(ID>=0 && ID<node_num){
                if(PartiTag[ID].first==-1){
                    PartiTag[ID].first=k;
                    PartiVertex[k].emplace_back(ID);
                }else{
                    cout<<"vertex already in one partition!"<<ID<<" "<<PartiTag[ID].first<<" "<<k<<endl;
                }
            }else{
                cout<<"Wrong vertex ID! "<<ID<<endl; exit(1);
            }

        }
    }
    //further check that each vertex is in one and only one partition
    for(int vid=0;vid<node_num;vid++){
        if(PartiTag[vid].first==-1){
            cout<<"vertex "<<vid<<" not within any partition"<<endl; exit(1);
        }
    }
    int nNum=0;
    for(int pid=0;pid<partiNum;++pid){
        nNum+=PartiVertex[pid].size();
    }
    if(nNum!=node_num){
        cout<<"Inconsistent node number! "<<nNum<<" "<<node_num<<endl; exit(1);
    }
    //record the vertex to PartiVertex in vertex order: from lower-rank vertex to higher-rank vertex

    ifstream IF(filename+"/subgraph_edge");
    if(!IF){
        cout<<"Cannot open file "<<"subgraph_edge"<<endl;
        exit(1);
    }

    int pnum1;
    IF>>pnum1;
    for(int k=0;k<pnum1;k++){
        int edgenum0,ID1,ID2,weight;
        IF>>edgenum0;
        for(int i=0;i<edgenum0;i++){
            IF>>ID1>>ID2>>weight;
//            if(flag_minus){
//                ID1 = ID1-1; ID2 = ID2-1;
//            }
            if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
                NeighborsParti[ID1].emplace_back(ID2,weight);
            }else{
                cout<<"Wrong for subgraph_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
            }


        }
    }

    vector<int> vecint;
    vecint.clear();
    NeighborSketch.assign(partiNum, vecint);
//    vector<set<int>> NeighborSketchS;
    NeighborSketchS.assign(partiNum, unordered_set<int>());

    BoundVertex.assign(partiNum,vector<vertex>());
    //read the cut edges
    ifstream IF2(filename+"/cut_edges");
    if(!IF2){
        cout<<"Cannot open file "<<"cut_edges"<<endl;
        exit(1);
    }

    int ednum,ID1,ID2,weight;
    int boundaryNum=0;
    int PID1, PID2;

    IF2>>ednum;
    for(int i=0;i<ednum;i++){
        IF2>>ID1>>ID2>>weight;
//        if(flag_minus){
//            ID1 = ID1-1; ID2 = ID2-1;
//        }

        if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
            PID1=PartiTag[ID1].first, PID2=PartiTag[ID2].first;
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                cout<<"two end points of cut edge are in the same partition"<<endl; exit(1);
            }
            if(!PartiTag[ID1].second){
                PartiTag[ID1].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID1].first].emplace_back(ID1);
            }
            if(!PartiTag[ID2].second){
                PartiTag[ID2].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID2].first].emplace_back(ID2);
            }

            NeighborsOverlay[ID1].insert({ID2,weight});
//            NeighborsOverlay[ID1].emplace_back(ID2,weight);

            if(NeighborSketchS[PID1].find(PID2)==NeighborSketchS[PID1].end()){//if not found, i.e., PID2 is not in the NeighborSketchS of PID1
                NeighborSketch[PID1].push_back(PID2);
                NeighborSketch[PID2].push_back(PID1);

                NeighborSketchS[PID1].insert(PID2);
                NeighborSketchS[PID2].insert(PID1);
            }
        }else{
            cout<<"Wrong for cut_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }
    }

    nodeNumOverlay=boundaryNum;


    /*for(int k=0;k<pnum;k++){
        cout<<k<<" "<<NeighborSketch[k].size()<<endl;
    }*/
}

void Graph::ReadPartitionGraph(){
    NeighborsParti.assign(node_num, vector<pair<vertex,int>>());
    NeighborsOverlay.assign(node_num,unordered_map<vertex,int>());
//    NeighborsOverlay.assign(node_num,vector<pair<vertex,int>>());
    PartiTag.assign(node_num, make_pair(-1,false));

    bool flag_minus = false;

    string filename=sourcePath+"/partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    ifstream IF1(filename+"/subgraph_vertex");
    if(!IF1){
        cout<<"Cannot open file "<<"subgraph_vertex"<<endl;
        exit(1);
    }

    int pnum2;
    IF1>>pnum2;
    if(algoParti == "NC"){
//        flag_minus = true;
        partiNum = pnum2;
    }else if(algoParti == "SC" || algoParti == "MT"){
//        flag_minus = true;
//        pnum2 = pnum;
    }
    cout<<"Partition number: "<<pnum2<<endl;

    PartiVertex.assign(partiNum,vector<vertex>());
    for(int k=0;k<pnum2;k++){
        int vernum,ID;
        IF1>>vernum;
        for(int i=0;i<vernum;i++){
            IF1>>ID;
//            if(flag_minus){
//                ID = ID-1;
//            }

            if(ID>=0 && ID<node_num){
                if(PartiTag[ID].first==-1){
                    PartiTag[ID].first=k;
                    PartiVertex[k].emplace_back(ID);
                }else{
                    cout<<"vertex already in one partition!"<<ID<<" "<<PartiTag[ID].first<<" "<<k<<endl;
                }
            }else{
                cout<<"Wrong vertex ID! "<<ID<<endl; exit(1);
            }

        }
    }
    //further check that each vertex is in one and only one partition
    for(int vid=0;vid<node_num;vid++){
        if(PartiTag[vid].first==-1){
            cout<<"vertex "<<vid<<" not within any partition"<<endl; exit(1);
        }
    }
    int nNum=0;
    for(int pid=0;pid<partiNum;++pid){
        nNum+=PartiVertex[pid].size();
    }
    if(nNum!=node_num){
        cout<<"Inconsistent node number! "<<nNum<<" "<<node_num<<endl; exit(1);
    }
    //record the vertex to PartiVertex in vertex order: from lower-rank vertex to higher-rank vertex

    ifstream IF(filename+"/subgraph_edge");
    if(!IF){
        cout<<"Cannot open file "<<"subgraph_edge"<<endl;
        exit(1);
    }

    int pnum1;
    IF>>pnum1;
    for(int k=0;k<pnum1;k++){
        int edgenum0,ID1,ID2,weight;
        IF>>edgenum0;
        for(int i=0;i<edgenum0;i++){
            IF>>ID1>>ID2>>weight;
//            if(flag_minus){
//                ID1 = ID1-1; ID2 = ID2-1;
//            }
            if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
                NeighborsParti[ID1].emplace_back(ID2,weight);
            }else{
                cout<<"Wrong for subgraph_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
            }


        }
    }

    BoundVertex.assign(partiNum,vector<vertex>());
    //read the cut edges
    ifstream IF2(filename+"/cut_edges");
    if(!IF2){
        cout<<"Cannot open file "<<"cut_edges"<<endl;
        exit(1);
    }

    int ednum,ID1,ID2,weight;
    int boundaryNum=0;
    int PID1, PID2;

    IF2>>ednum;
    for(int i=0;i<ednum;i++){
        IF2>>ID1>>ID2>>weight;
//        if(flag_minus){
//            ID1 = ID1-1; ID2 = ID2-1;
//        }

        if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
            PID1=PartiTag[ID1].first, PID2=PartiTag[ID2].first;
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                cout<<"two end points of cut edge are in the same partition"<<endl; exit(1);
            }
            if(!PartiTag[ID1].second){
                PartiTag[ID1].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID1].first].emplace_back(ID1);
            }
            if(!PartiTag[ID2].second){
                PartiTag[ID2].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID2].first].emplace_back(ID2);
            }

            NeighborsOverlay[ID1].insert({ID2,weight});
//            NeighborsOverlay[ID1].emplace_back(ID2,weight);

        }else{
            cout<<"Wrong for cut_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }
    }

    nodeNumOverlay=boundaryNum;


}

void Graph::OverlayOrderingBuild(){

    int lastParti = -1;
    map<int,pair<int,int>> m;
    E.assign(partiNum,m);
    for(int i=0;i<NeighborSketch.size();i++){
        for(auto it=NeighborSketch[i].begin();it!=NeighborSketch[i].end();++it){
            E[i].insert(make_pair(*it,make_pair(1,1)));
        }
    }
    _DD_.assign(partiNum,0);
    DD.assign(partiNum,0);
    _DD2_.assign(partiNum,0);

    set<DegComp> Deg;
    set<DegComp2> Deg2;
    int ID,degree;
    for(ID=0;ID<partiNum;ID++){
        degree=NeighborSketch[ID].size();
        if(degree!=0){
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
        }else{
            cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
        }
    }

    vector<bool> exist; exist.assign(partiNum,true);
    vector<bool> change; change.assign(partiNum,false);

    vector<set<int>> NeighborCon(partiNum,set<int>());
    vector<int> neix;
    neix.clear();

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    int x;

    set<int> pSet;
    while(!Deg.empty()){

        x=(*Deg.begin()).x;
        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }
        Deg.erase(Deg.begin());



        if(lastParti!=-1 && NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if x is lastParti's neighbor
            neix.clear();
            while(NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if found
                _DD2_[x]++;
                neix.emplace_back(x);

                if(Deg.empty()){
                    break;
                }else{
                    x=Deg.begin()->x;
                    Deg.erase(Deg.begin());
                }
            }

            if(NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if found, i.e., x is the neighbor of lastParti
                if(neix.size()>1){//if neix has more than one element
                    if(Deg.empty()){
                        Deg2.clear();
                        for(int i=0;i<neix.size();++i){
                            Deg2.insert(DegComp2(neix[i]));
                        }
                        x=Deg2.begin()->x;
                        Deg2.erase(Deg2.begin());
                        if(!Deg2.empty()){
                            for(auto it=Deg2.begin();it!=Deg2.end();++it){
                                ID=it->x;
                                Deg.insert(DegComp(ID));
                            }
                            Deg2.clear();
                        }
                    }else{
                        cout<<"Wrong! "<<endl; exit(1);
                    }
                }
            }//if not the neighbor
            else{
                if(!neix.empty()){
                    for(int i=0;i<neix.size();++i){
                        Deg.insert(DegComp(neix[i]));
                    }
                }
            }
        }

//        cout<<x<<" "<<Deg.size()<<endl;
        vNodeOrderOverlay.emplace_back(x);
        pSet.insert(x);
        lastParti = x;
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.emplace_back(*it);
                NeighborCon[x].insert(it->first);
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderOverlay.size() != partiNum || pSet.size()!=partiNum){
        cout<<"Inconsistent size for sketch graph! "<<vNodeOrderOverlay.size()<<" "<<pSet.size() <<" "<< partiNum<<endl; exit(1);
    }

//    exit(0);
}
//boundary-first heuristic ordering, minimize overlay tree height
void Graph::OverlayOrderingBoundaryFirstMDEGreedy(){
    for(int i=0;i<NeighborsOverlay.size();i++){
        if(!NeighborsOverlay[i].empty()){
            for(auto it=NeighborsOverlay[i].begin();it!=NeighborsOverlay[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
        }
    }

    set<DegHeightComp> Deg;
    int ID,degree;
    for(int pid=0;pid<partiNum;++pid){
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID = BoundVertex[pid][i];
            degree=E[ID].size();
            if(degree!=0){
                _DD_[ID]=degree;
                DD[ID]=degree;
                Deg.insert(DegHeightComp(ID));
            }else{
                cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
            }
        }
    }


//    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    int maxHeight=0;
    cout<<"Overlay vertex number: "<<Deg.size()<<endl;
    Timer tt;
    tt.start();
    int stepShow=Deg.size()/10;
    while(!Deg.empty()){
//        if(count%stepShow==0){
//            tt.stop();
//            cout<<"count "<<count<<" , treewidth "<<Twidth<<" "<< tt.GetRuntime()<<endl;
//            tt.start();
//        }


        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(Deg.begin());
                _DD_[x]=DD[x];
                Deg.insert(DegHeightComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }
        count+=1;
//        vNodeOrderParti[pid].insert({order_i,x});
        vNodeOrderOverlay.push_back(x);
        order_i++;
        Deg.erase(Deg.begin());

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            Neigh.emplace_back(*it);
            if(_DD2_[it->first]<_DD2_[x]+1){//obtain the estimated tree height
                _DD2_[it->first]=_DD2_[x]+1;
//                change[it->first]=true;
            }
        }
        if(maxHeight<_DD2_[x]){
            maxHeight=_DD2_[x];
        }
        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    cout<<"Max tree height: "<<maxHeight<<" ; treewidth: "<<Twidth<<endl;
//    if(vNodeOrderOverlay.size() != partiNum || pSet.size()!=partiNum){
//        cout<<"Inconsistent size for sketch graph! "<<vNodeOrderOverlay.size()<<" "<<pSet.size() <<" "<< partiNum<<endl; exit(1);
//    }

//    exit(0);
}
//boundary-first heuristic ordering, simple MDE
/*void Graph::OverlayOrderingBoundaryFirstMDE(){
    for(int i=0;i<NeighborsOverlay.size();i++){
        if(!NeighborsOverlay[i].empty()){
            for(auto it=NeighborsOverlay[i].begin();it!=NeighborsOverlay[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
        }
    }

    set<DegComp> Deg;
    int ID,degree;
    for(int pid=0;pid<partiNum;++pid){
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID = BoundVertex[pid][i];
            degree=E[ID].size();
            if(degree!=0){
                _DD_[ID]=degree;
                DD[ID]=degree;
                Deg.insert(DegComp(ID));
            }else{
                cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
            }
        }
    }


//    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    int maxHeight=0;
    cout<<"Overlay vertex number: "<<Deg.size()<<endl;
    Timer tt;
    tt.start();
    while(!Deg.empty()){
        if(count%500==0){
            tt.stop();
            cout<<"count "<<count<<" , treewidth "<<Twidth<<" "<< tt.GetRuntime()<<endl;
            tt.start();
        }

        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(Deg.begin());
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }
        count+=1;
//        vNodeOrderParti[pid].insert({order_i,x});
        vNodeOrderOverlay.push_back(x);
        order_i++;
        Deg.erase(Deg.begin());

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            Neigh.emplace_back(*it);
            if(_DD2_[it->first]<_DD2_[x]+1){//obtain the estimated tree height
                _DD2_[it->first]=_DD2_[x]+1;
//                change[it->first]=true;
            }
        }
        if(maxHeight<_DD2_[x]){
            maxHeight=_DD2_[x];
        }
        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    cout<<"Max tree height: "<<maxHeight<<" ; treewidth: "<<Twidth<<endl;
//    if(vNodeOrderOverlay.size() != partiNum || pSet.size()!=partiNum){
//        cout<<"Inconsistent size for sketch graph! "<<vNodeOrderOverlay.size()<<" "<<pSet.size() <<" "<< partiNum<<endl; exit(1);
//    }

//    exit(0);
}*/
//based on SketchOrder function, to guarantee the order of neighboring vertices are not contiguous
void Graph::OverlayOrderingBuildBoundaryFirst(int nodenum, vector<vector<int>> Neighbor){
    map<int,pair<int,int>> m;
    E.assign(partiNum,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j],make_pair(1,1)));
    }

    _DD_.assign(partiNum,0);
    DD.assign(partiNum,0);

    set<DegComp> Deg;
    int degree;
    for(int i=0;i<Neighbor.size();i++){
        degree=Neighbor[i].size();
        if(degree!=0){
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
        }
    }

    vector<bool> exist; exist.assign(partiNum,true);
    vector<bool> change; change.assign(partiNum,false);

    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(partiNum,vect); //NeighborCon.clear();
    //SCconNodes.clear();

    //cout<<"Begin to contract"<<endl;
    int count=0;
    int x;
    int lastx;
    vector<int> neix;
    neix.clear();
    while(!Deg.empty()){
        //if(count%10==0) cout<<"count "<<count<<endl;
        count+=1;

        while(true){
            if(Deg.empty()){//in case that all the remaining vertices are the neighbors of current vertex x
                x=neix[0];
                for(int j=1;j<neix.size();j++){
                    Deg.insert(DegComp(neix[j]));
//					cout<<"insert back/// "<<neix[j]<<endl;
                }
                neix.clear();
                break;///
            }
            else
                x=(*Deg.begin()).x;

            while(true){
                if(change[x]){
                    Deg.erase(DegComp(x));
                    _DD_[x]=DD[x];
                    Deg.insert(DegComp(x));
                    change[x]=false;
                    x=(*Deg.begin()).x;
                }else
                    break;
            }

            if(count==1){
                lastx=x;
                break;
            }else if(NeighborSketchS[lastx].find(x)==NeighborSketchS[lastx].end()){//if not found
                lastx=x;
                break;
            }else{
                Deg.erase(DegComp(x));
                neix.push_back(x);
//				cout<<"erase "<<x<<endl;
            }

//            if(Deg.empty())////
//                break;
        }

        if(neix.size()!=0){
            for(int j=0;j<neix.size();j++){
                Deg.insert(DegComp(neix[j]));
//				cout<<"insert back "<<neix[j]<<endl;
            }
        }
        neix.clear();

        vNodeOrderOverlay.push_back(x);
        Deg.erase(DegComp(x));
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();i++){
            for(int j=i+1;j<Neigh.size();j++){
                insertEOrderGenerate(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                change[Neigh[i].first]=true;
                change[Neigh[j].first]=true;
            }
        }
    }

    /*NodeOrderSketch.assign(nodenum,-1);
    for(int k=0;k<vNodeOrderSketch.size();k++){
        NodeOrderSketch[vNodeOrderSketch[k]]=k;
        cout<<"order "<<k<<", vertex "<<vNodeOrderSketch[k]<<", degree "<<NeighborSketch[vNodeOrderSketch[k]].size()<<endl;
    }*/
    //cout<<"Finish Contract"<<endl;
}
//MDE-based vertex ordering for partition
void Graph::PartitionOrderingBuildMDE(bool ifParallel){
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsParti.size();i++){
        for(auto it=NeighborsParti[i].begin();it!=NeighborsParti[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    _DD_.assign(node_num,0);
    DD.assign(node_num,0);
//    vNodeOrderParti.assign(partiNum,map<int,int>());
    vNodeOrderParti.assign(partiNum,vector<int>());

    if(ifParallel){
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrderingV, this, boost::ref(processID[j])));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrdering, this, j));
            }
            thread.join_all();
        }
    }
    else{
        for(int i=0;i<partiNum;++i){
            PartitionOrdering(i);
        }
    }

}

//Boundary-first MDE-based vertex ordering for partition
void Graph::PartitionOrderingBoundaryFirstMDE(bool ifParallel, vector<bool>& change){
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsParti.size();i++){
        for(auto it=NeighborsParti[i].begin();it!=NeighborsParti[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    _DD_.assign(node_num,0);
    _DD2_.assign(node_num,1);//estimated tree height
    DD.assign(node_num,0);
//    vNodeOrderParti.assign(partiNum,map<int,int>());
    vNodeOrderParti.assign(partiNum,vector<int>());
    heightMaxs.assign(partiNum,1);
    if(ifParallel){
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrderingBoundaryFirstV, this, boost::ref(processID[j]), boost::ref(change) ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrderingBoundaryFirst, this, j, boost::ref(change) ));
            }
            thread.join_all();
        }
    }
    else{
        for(int i=0;i<partiNum;++i){
            PartitionOrderingBoundaryFirst(i, change);
        }
    }

    double aveHeight=0.0;
    int maxHeight=0;
    for(int i=0;i<heightMaxs.size();++i){
        aveHeight+=heightMaxs[i];
        if(maxHeight<heightMaxs[i]){
            maxHeight=heightMaxs[i];
        }
    }
    cout<<"Maximal in-partition tree height: "<<maxHeight<<" ; average in-partition tree height: "<<aveHeight/partiNum<<" ; boundary-clique partition number: "<<cliqueBoundaryPNum<<endl;
}

void Graph::PartitionOrderingV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        PartitionOrdering(p[i]);
    }
}

void Graph::PartitionOrdering(int pid){

    set<DegComp> Deg;
    int ID,degree;
    for(int i=0;i<PartiVertex[pid].size();i++){
        ID = PartiVertex[pid][i];
        degree=NeighborsParti[ID].size();
        if(degree!=0){
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
        }else{
            cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
        }
    }

    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    while(!Deg.empty()){
//        if(count%10000==0)
//            cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
        count+=1;
        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }

//        vNodeOrderParti[pid].insert({order_i,x});
        vNodeOrderParti[pid].push_back(x);
        order_i++;
        Deg.erase(Deg.begin());
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.emplace_back(*it);
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderParti[pid].size() != PartiVertex[pid].size()){
        cout<<"Inconsistent size! "<< pid <<" "<<vNodeOrderParti[pid].size() <<" "<< PartiVertex[pid].size()<<endl; exit(1);
    }
}

void Graph::PartitionOrderingBoundaryFirstV(vector<int>& p, vector<bool>& change){
    for(int i=0;i<p.size();++i){
        PartitionOrderingBoundaryFirst(p[i], change);
    }
}

void Graph::PartitionOrderingBoundaryFirst(int pid, vector<bool>& change){

    set<DegComp> Deg;
    int ID,degree;
    for(int i=0;i<PartiVertex[pid].size();i++){
        ID = PartiVertex[pid][i];
        if(!PartiTag[ID].second) {//if not boundary
            degree = NeighborsParti[ID].size();
            if(degree!=0){
                _DD_[ID]=degree;
                DD[ID]=degree;
                Deg.insert(DegComp(ID));
            }else{
                cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
            }
        }else{//skip boundary vertex
            degree = NeighborsParti[ID].size();
            if(degree!=0){
                _DD_[ID]=degree;
                DD[ID]=degree;
            }else{
                cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
            }
        }

    }

//    vector<bool> exist; exist.assign(node_num,true);
//    vector<bool> change; change.assign(node_num,false);

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    int heightTemp=0;
    while(!Deg.empty()){
//        if(count%10000==0)
//            cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
        count+=1;
        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(Deg.begin());
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }

//        vNodeOrderParti[pid].insert({order_i,x});
        vNodeOrderParti[pid].push_back(x);
        order_i++;
        Deg.erase(Deg.begin());

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            Neigh.emplace_back(*it);
            if(_DD2_[it->first]<_DD2_[x]+1){//obtain the estimated tree height
                _DD2_[it->first]=_DD2_[x]+1;
            }
            if(!PartiTag[it->first].second ){//if not boundary vertex
                if(heightTemp<_DD2_[it->first]){
//                    cout<<pid<<" : "<<it->first<<" , "<<x<<" : "<<heightTemp<<" "<<treeHeight2[it->first]<<endl;
                    heightTemp=_DD2_[it->first];
                }
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderParti[pid].size()+BoundVertex[pid].size() != PartiVertex[pid].size()){
        cout<<"Inconsistent size! "<< pid <<" "<<vNodeOrderParti[pid].size() <<" "<< PartiVertex[pid].size()<<endl; exit(1);
    }
    //obtain the tree height
//    vector<int> treeHeight(node_num,INF);
//    vector<int> rankTemp(node_num,INF);
//    ID=vNodeOrderParti[pid][vNodeOrderParti[pid].size()-1];
//    rankTemp[ID]=1;
//    treeHeight[ID]=1;
//    heightMaxs[pid]=1;
//    for(int i=vNodeOrderParti[pid].size()-2;i>-1;--i){
//        ID=vNodeOrderParti[pid][i];
//        rankTemp[ID]=vNodeOrderParti[pid].size()-i;
//
//        if(!E[ID].empty()){
//            bool ifRoot=true;
//            int pa=-1;
//            for(auto it=E[ID].begin();it!=E[ID].end();++it){
//                if(!PartiTag[it->first].second){//not boundary vertex
//                    if(pa==-1){
//                        pa=it->first;
//                    }else{
//                        if(rankTemp[it->first]>rankTemp[pa]){//pick the vertex with the lowest order as the parent
//                            pa=it->first;
//                        }
//                    }
//                    ifRoot=false;
//                }
//            }
//            if(ifRoot){
////                cout<<"New root by all boundary neighbor: "<<ID<<endl;
//                treeHeight[ID]=1;
//            }else{
//                treeHeight[ID]=treeHeight[pa]+1;
//            }
//
//        }
//        else{
//            cout<<"New root by empty neighbor: "<<ID<<endl;
//            treeHeight[ID]=1;
//        }
////        cout<<ID<<" "<<i<<" "<<rankTemp[ID]<<" "<<treeHeight[ID]<<endl;
//        if(heightMaxs[pid]<treeHeight[ID]){
////            cout<<pid<<" : "<<heightMaxs[pid]<<" "<<treeHeight[ID]<<endl;
//            heightMaxs[pid]=treeHeight[ID];
//        }
//    }
//    if(heightMaxs[pid]!=heightTemp){
//        cout<<pid<<" Inconsistent tree height: "<<heightTemp+1<<" "<<heightMaxs[pid]<<endl;
//    }
    heightMaxs[pid]=heightTemp;
    int cliqueNum=0;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID=BoundVertex[pid][i];
        if(E[ID].size()==BoundVertex[pid].size()-1){
            cliqueNum++;
        }
//        else{
//            cout<<pid<<": "<<ID<<" "<<E[ID].size()<<" "<<BoundVertex[pid].size()<<endl;
//        }
    }
    if(cliqueNum==BoundVertex[pid].size()){
//        cout<<"Partition "<<pid<<" is not a clique. "<<cliqueNum<<" "<<BoundVertex[pid].size()<<endl;
        sm->wait();
        cliqueBoundaryPNum++;
        sm->notify();
    }

}



/// Functions for MDE contraction
void Graph::deleteEOrderGenerate(int u,int v){//delete u from v's neighbor
    // Unnecessary erase
//    if(E[u].find(v)!=E[u].end()){
//        E[u].erase(E[u].find(v));
//        DD[u]--;
//    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}

void Graph::insertEOrderGenerate(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//        DD2[u]++;
    }
    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        DD[v]++;
//        DD2[u]++;
    }
}


/// Query Processing
//function for correctness check
void Graph::CorrectnessCheck(int runtimes){
    Timer tt;
    double runT=0;
    srand (time(NULL));
    int s, t, d1, d2, d3;
//    runtimes = 1;
//    algoQuery=PCH_No;
//    algoQuery=PH2H_Post;
    cout<<"Correctness check ("<<runtimes<<" rounds) ... AlgoQuery: "<<algoQuery<<". ";
    int sameNum=0, crossNum=0;
    innerBenefitNum=0;
    outerBenefitNum=0;
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
//        int pid=rand()%partiNum;
//        s=PartiVertex[pid][rand()%PartiVertex[pid].size()];
//        t=PartiVertex[pid][rand()%PartiVertex[pid].size()];

//        s=859796,t=546223;//VPL, DCH, cross
//        s=706838,t=520171;//VPL, post
//        s=935728,t=697993;//VPL, cross
//        s=427834,t=697993;//VPL, cross
//        s=457100,t=427834;//VPL, cross


//        s=400912,t=213579;//VPL
//        s=213042,t=821345;//VPL, Cross
//        s=704986,t=691421;//VPL

//        s=765745,t=629253;//VPL
//        s=769156,t=632036;//VPL
//        s=764513,t=630166;//VPL
//        s=629367,t=630166;//VPL
//        s=632597,t=629367;//VPL
//        s=633083,t=629367;//VPL
//        cout<<"Query "<<i<<": "<<s<<" "<<t<<endl;

//        if(runtimes == 1){
//            cout<<"s: "<<s<<" ; t: "<<t<<endl;
//        }

        if(algoChoice==1){
            tt.start();
            d2=QueryNP(s,t);
            tt.stop();
        }else if(algoChoice==2){
            tt.start();
            d2=Query(s,t);
            tt.stop();
        }else if(algoChoice==3){
            tt.start();
            d2=QueryPMHL(s,t);
            tt.stop();
            if(PartiTag[s].first==PartiTag[t].first){
                sameNum++;
            }else{
                crossNum++;
            }
        }
        else if(algoChoice==4){
            tt.start();
            d2=QueryPMHLOpt(s,t);
            tt.stop();
        }
        else if(algoChoice==5){
            tt.start();
            d2=QueryPostMHL(s,t);
            tt.stop();
        }
        else if(algoChoice==6){
            tt.start();
//            d2=QueryVPL(s,t);
            d2=QueryVPL_HASP(s,t);
//            d2=QueryVPL_HASPDebug(s,t);
            tt.stop();
            if(PartiTag[s].first==PartiTag[t].first){
                sameNum++;
            }else{
                crossNum++;
            }
        }
        else if(algoChoice==0){
            tt.start();
            d2= Astar(s,t,Neighbor);
            tt.stop();
        }
        runT+=tt.GetRuntime();

        d1=Dijkstra(s,t,Neighbor);
//        cout<<"Algorithm "<<algoQuery<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTag[s].first<<" "<<PartiTag[t].first<<"; Boundary Tag: "<<PartiTag[s].second<<" "<<PartiTag[t].second<<endl;

//        cout<<i<<": "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl;
        if(d1!=d2){
            if(algoChoice==1){
                cout<<"Correct Test. InCorrect! Algorithm "<<algoQuery<<", "<<i<<": "<<s<<" "<<t<<": "<<d2<<" "<<d1<<endl;
            }else if(algoChoice==2 || algoChoice==3 || algoChoice==4){
                cout<<"Correct Test. InCorrect! Algorithm "<<algoQuery<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTag[s].first<<" "<<PartiTag[t].first<<"; Boundary Tag: "<<PartiTag[s].second<<" "<<PartiTag[t].second<<endl;
                QueryDebug(s,t);
            }else if(algoChoice==5){
                cout<<"Correct Test. InCorrect! Algorithm "<<algoQuery<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTags[s].first<<" "<<PartiTags[t].first<<endl;
                QueryPostMHLDebug(s,t);
            }else if(algoChoice==6){
                cout<<"Correct Test. InCorrect! Algorithm "<<algoQuery<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTag[s].first<<" "<<PartiTag[t].first<<"; Boundary Tag: "<<PartiTag[s].second<<" "<<PartiTag[t].second<<endl;
//                QueryVPLDebug(s,t);
                QueryVPL_HASPDebug(s,t);
            }else if(algoChoice==0){
                cout<<"InCorrect! "<<i<<": "<<s<<" "<<t<<": "<<d2<<" "<<d1<<endl;
            }
            exit(1);
        }
    }
    if(algoChoice==3 || algoChoice==6){
        cout<<"Same-partition query number: "<<sameNum<<" ; Cross-partition query number: "<<crossNum<<endl;
    }
    if(algoChoice==6){
        cout<<"Inner-unaffected query number: "<<innerBenefitNum<<" ; Outer-unaffected query number: "<<outerBenefitNum<<endl;
    }
    cout<<" Average Query Time: "<<1000*runT/runtimes<<" ms."<<endl;
}

//function for Query processing, debug version
int Graph::QueryDebug(int ID1, int ID2){
    int dis=INF;

    if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in core
        cout<<"Core-Core"<<endl;
//        dis=QueryCoreDebug(ID1, ID2);

    }else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
        cout<<"Core-Parti"<<endl;
        dis=QueryPartiCoreDebug(ID2, ID1);

    }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
        cout<<"Parti-Core"<<endl;
        dis=QueryPartiCoreDebug(ID1, ID2);

    }else if(!PartiTag[ID1].second && !PartiTag[ID2].second){//both in partition

        if(PartiTag[ID1].first != PartiTag[ID2].first){//Case 3: in different peripheries
            cout<<"Parti-Parti"<<endl;
            int d=INF;
            int b1,b2,d1,d2;//final results
            int pid1=PartiTag[ID1].first;
            int pid2=PartiTag[ID2].first;

            vector<int> B1=BoundVertex[pid1];
            vector<int> B2=BoundVertex[pid2];

            map<int,int> m1,m2;
            m1.clear();
            m2.clear();
            int bID1, bID2, tempdis;
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];

//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
                m1.insert(make_pair(bID1, QueryH2HPartition(ID1,bID1,pid1)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
                m2.insert(make_pair(bID2, QueryH2HPartition(ID2,bID2,pid2)));
            }

            for(int k=0;k<B1.size();k++){
                bID1=B1[k];

                if(m1[bID1]>d)
                    continue;

                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];

                    if(m2[bID2]>d)
                        continue;

                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d){
                        d=tempdis;
                        b1=bID1; b2=bID2; d1=m1[bID1]; d2=m2[bID2];
//                        cout<<b1<<" "<<b2<<" "<<d<<endl;
                    }
                }
            }
            dis=d;
            int d_12=QueryCore(b1,b2), dDijk_s=Dijkstra(ID1,b1,Neighbor), dDijk_12=Dijkstra(b1,b2,Neighbor), dDijk_t=Dijkstra(b2,ID2,Neighbor);
            cout<<ID1<<" "<<b1<<"("<<NodeOrder[b1]<<","<<PartiTag[b1].first<<") "<<b2<<"("<<NodeOrder[b2]<<","<<PartiTag[b2].first<<") "<<ID2<<" : "<<d1<<" "<<d_12<<" "<<d2<<" ; "<<dDijk_s<<" "<<dDijk_12<<"("<<DijkstraCore(b1,b2)<<") "<<dDijk_t<<endl;

            QueryCoreDebug(b1,b2);
            QueryPartiPartiExtLCADebug(ID1,ID2);

//                if(d1!=dDijk_s){
//                    DijkstraPath(ID1,b1);
//                }
//                if(d_12!=dDijk_12){
//                    DijkstraPath(b1,b2);
//                }
//                if(d2!=dDijk_t){
//                    DijkstraPath(b2,ID2);
//                }

        }
        else{//Case 4: in the same periphery
            cout<<"Same-Parti"<<endl;
//                dis= QuerySameParti(ID1,ID2);
            int d=INF;
            int b1,b2,df1,df2;
            int pid1=PartiTag[ID1].first;
            int pid2=PartiTag[ID2].first;

            int temp_dis = QueryH2HPartition(ID1, ID2, pid1);/// d2 may be wrong sometimes
            if (temp_dis < d){
                d = temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
                b1=b2=-1;
                df1=df2=-1;
            }

            vector<int> B = BoundVertex[pid1];
            map<int, int> m1, m2;
            m1.clear();
            m2.clear();
            vector<int> B1, B2;
            B1.clear();
            B2.clear();
            int bID, d1, d2;
            for (int i = 0; i < B.size(); i++) {
                bID = B[i];

                d1 = QueryH2HPartition(ID1,bID,pid1);
                d2 = QueryH2HPartition(ID2,bID,pid1);

                if (d1 < d) {
                    B1.push_back(bID);
                    m1.insert(make_pair(bID, d1));
                }
                if (d2 < d) {
                    B2.push_back(bID);
                    m2.insert(make_pair(bID, d2));
                }
            }

            int bID1, bID2, tempdis;
            if (!B1.empty() && !B2.empty()) {
                for (int k = 0; k < B1.size(); k++) {
                    bID1 = B1[k];
                    if (m1[bID1] > d)
                        continue;
                    for (int z = 0; z < B2.size(); z++) {
                        bID2 = B2[z];
                        if (m2[bID2] > d)
                            continue;
                        tempdis = m1[bID1] + QueryCore(bID1, bID2) + m2[bID2];
                        if (tempdis < d){
                            d = tempdis;
                            b1=bID1;b2=bID2;
                            df1=m1[bID1];df2=m2[bID2];
                        }
                    }
                }
            }

            if(b1!=-1){
                cout<<"d4: "<<ID1<<" "<<b1<<" "<<b2<<" "<<ID2<<" : "<<df1<<" "<<QueryCore(b1,b2)<<" "<<df2<<" ; "<<Dijkstra(ID1,b1,Neighbor)<<" "<<Dijkstra(b1,b2,Neighbor)<<" "<<Dijkstra(b2,ID2,Neighbor)<<endl;
            }else{
                int dDijk2 = Dijkstra(ID1,ID2,Neighbor);
                cout<<"d2: "<<d<<"; "<<dDijk2<<endl;
                if(d!=dDijk2){
//                        DijkstraPath(ID1,ID2);
                }
            }

            dis = d;

        }

    }
    return dis;
}

//function for core index correctness check
void Graph::CorrectnessCheckCore(int runtimes){
    srand (time(NULL));
    int s, t, d1, d2, d3;
    vector<int> coreVertex;
    for(int i=0;i<node_num;++i){
        if(PartiTag[i].second){
            coreVertex.emplace_back(i);
        }
    }
    int corenum=coreVertex.size();
    cout<<"Core graph correctness check ("<<runtimes<<" rounds)..."<<endl;
    for(int i=0;i<runtimes;i++){
        s=coreVertex[rand()%corenum];
        t=coreVertex[rand()%corenum];
        if(PartiTag[s].second && PartiTag[t].second){//for core vertex
//            d1=QueryCore(s,t);
            d1=DijkstraCore(s,t);
            d2=Dijkstra(s,t,Neighbor);

            if(d1!=d2){
                cout<<"InCorrect! "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<"): "<<d1<<" "<<d2<<endl;
//				DijkstraPath(s,t);
//				DijkstraCorePath(s,t);
                exit(1);
            }
        }else
            i--;
    }
}

//function for efficiency test
void Graph::EffiCheck(string filename,int runtimes){
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Query file: "<<filename<<endl;
    int num, ID1, ID2;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.push_back(make_pair(ID1, ID2));
    }
    if(runtimes > num){
        runtimes = num;
    }
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;

    double runT=0;
    int d1,d2;
    bool ifDebug=false;
//    ifDebug=true;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }
    clock_t start = clock();
    vector<int> results(runtimes,-1);
    int crossNum=0, sameNum=0;
    for(int i=0;i<runtimes;i++){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
//        if(PartiTag[ID1].first!=PartiTag[ID2].first){
//            cout<<"Different Partition: "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<endl;
//        }
        if(algoChoice==1){
            tt.start();
            d1=QueryNP(ID1,ID2);
            tt.stop();
        }else if(algoChoice==2){
            tt.start();
            d1=Query(ID1,ID2);
            tt.stop();
        }else if(algoChoice==3){
            tt.start();
            d1=QueryPMHL(ID1,ID2);
            tt.stop();
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                sameNum++;
            }else{
                crossNum++;
            }
        }else if(algoChoice==4){
            tt.start();
            d1=QueryPMHLOpt(ID1,ID2);
            tt.stop();
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                sameNum++;
            }else{
                crossNum++;
            }
        }else if(algoChoice==5){
            tt.start();
            d1=QueryPostMHL(ID1,ID2);
            tt.stop();
        }else if(algoChoice==6){
            tt.start();
            d1=QueryVPL(ID1,ID2);
            tt.stop();
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                sameNum++;
            }else{
                crossNum++;
            }
        }else if(algoChoice==0){
            tt.start();
            d1=Astar(ID1,ID2,Neighbor);
//            d1=Dijkstra(ID1,ID2,Neighbor);
            tt.stop();
        }

        runT+=tt.GetRuntime();
        results[i]=d1;

        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                if(algoChoice==5){
                    cout<<"Incorrect! "<<i<<": "<<ID1<<"("<<PartiTags[ID1].first<<") "<<ID2<<"("<<PartiTags[ID2].first<<"): "<<d1<<" "<<d2<<endl;
                    QueryPostMHLDebug(ID1,ID2);
                    exit(1);
                }else if(algoChoice>=1 && algoChoice<=4){
                    cout<<"Incorrect! "<<i<<": "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<"): "<<d1<<" "<<d2<<endl; exit(1);
                }

            }
        }

    }
    if(algoChoice==3 || algoChoice==6){
        cout<<"Same-partition query number: "<<sameNum<<" ; Cross-partition query number: "<<crossNum<<endl;
    }


    cout<<"Average Query Time: "<<(double)runT*1000/runtimes<<" ms."<<endl;
}

//function for efficiency test
void Graph::EffiCheckStages(vector<pair<int,int>> & ODpair, int runtimes, int intervalT, unsigned long long & throughputNum, vector<double>& stageUpdateT, vector<double>& stageQueryT){
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;

    double runT=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
    vector<double> qTime1, qTime2, qTime3, qTime4, qTime5;
    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    double updateT=0;//overall update time
    unsigned long long qt1=0,qt2=0,qt3=0,qt4=0,qt5=0,qtAll=0;
    int d1,d2;
    bool ifDebug=false;
//    ifDebug=true;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }
    clock_t start = clock();
    vector<int> results(runtimes,-1);

    if(algoChoice==1){
        runT1=EffiMHLStage(ODpair,1000,Dijk,qTime1);
        runT2=EffiMHLStage(ODpair,runtimes/10,CH,qTime2);
        runT3=EffiMHLStage(ODpair,runtimes,H2H,qTime3);
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT){
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: CH
            updateT+=stageDurations[CH];
            if(updateT<intervalT){
                dt2=stageDurations[CH]; qt2=dt2/runT2;
                //Stage 3: CH
                stageDurations[H2H]=intervalT-updateT;
                dt3=stageDurations[H2H]; qt3=dt3/runT3;
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }else{
            dt1=intervalT; qt1=dt1/runT1;
        }


        qtAll=qt1+qt2+qt3;
        cout<<"Stage 1 (BiDijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (CH): Duration: "<<dt2<<" ("<<stageDurations[CH]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (H2H): Duration: "<<dt3<<" ("<<stageDurations[H2H]<<") s; Throughput number: "<<qt3<<" ; query time: "<<1000 * runT3 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[CH]+=stageDurations[CH];
        stageUpdateT[H2H]+=stageDurations[H2H];
        stageQueryT[Dijk]+=runT1;
        stageQueryT[CH]+=runT2;
        stageQueryT[H2H]+=runT3;
    }
    else if(algoChoice==3){
        runT1=EffiPMHLStage(ODpair,1000,Dijk,qTime1);
        runT2=EffiPMHLStage(ODpair,runtimes/10,PCH_No,qTime2);
        runT3=EffiPMHLStage(ODpair,runtimes,PH2H_No,qTime3);
        runT4=EffiPMHLStage(ODpair,runtimes,PH2H_Post,qTime4);
        runT5=EffiPMHLStage(ODpair,runtimes,PH2H_Cross,qTime5);
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT){
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT){
                dt2=stageDurations[PCH_No]; qt2=dt2/runT2;
                //Stage 3: No-boundary PMHL
                updateT+=stageDurations[PH2H_No];
                if(updateT<intervalT){
                    dt3=stageDurations[PH2H_No]; qt3=dt3/runT3;
                    //Stage 4: Post-boundary PMHL
                    updateT+=stageDurations[PH2H_Post];
                    if(updateT<intervalT){
                        dt4=stageDurations[PH2H_Post]; qt4=dt4/runT4;
                        //Stage 5: Cross-boundary PMHL
                        stageDurations[PH2H_Cross]=intervalT-updateT;
                        dt5=stageDurations[PH2H_Cross]; qt5=dt5/runT5;
                    }else{
                        dt4=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]-stageDurations[PH2H_No]; qt4=dt4/runT4;
                    }
                }else{
                    dt3=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]; qt3=dt3/runT3;
                }
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }
        else{
            dt1=intervalT; qt1=dt1/runT1;
        }

        qtAll=qt1+qt2+qt3+qt4+qt5;
        cout<<"Stage 1 (Dijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (PCH): Duration: "<<dt2<<" ("<<stageDurations[PCH_No]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-No): Duration: "<<dt3<<" ("<<stageDurations[PH2H_No]<<") s; Throughput number: "<<qt3<<" ; query time: "<<1000 * runT3 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Post): Duration: "<<dt4<<" ("<<stageDurations[PH2H_Post]<<") s; Throughput number: "<<qt4<<" ; query time: "<<1000 * runT4 << " ms."<<endl;
        cout<<"Stage 5 (PH2H-Extend): Duration: "<<dt5<<" ("<<stageDurations[PH2H_Cross]<<") s; Throughput number: "<<qt5<<" ; query time: "<<1000 * runT5 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
        stageUpdateT[PH2H_No]+=stageDurations[PH2H_No];
        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
        stageQueryT[Dijk]+=runT1;
        stageQueryT[PCH_No]+=runT2;
        stageQueryT[PH2H_No]+=runT3;
        stageQueryT[PH2H_Post]+=runT4;
        stageQueryT[PH2H_Cross]+=runT5;
    }
    else if(algoChoice==5){
        runT1 = EffiPostMHLStage(ODpair,1000,Dijk,qTime1);
        runT2 = EffiPostMHLStage(ODpair, runtimes / 10, PCH_No, qTime2);
        runT4 = EffiPostMHLStage(ODpair, runtimes, PH2H_Post, qTime4);
        runT5 = EffiPostMHLStage(ODpair, runtimes, PH2H_Cross, qTime5);
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT) {
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT) {
                dt2=stageDurations[PCH_No]; qt2=dt2/runT2;
                //Stage 4: Post-boundary PMHL
                updateT+=stageDurations[PH2H_Post];
                if(updateT<intervalT) {
                    dt4=stageDurations[PH2H_Post]; qt4=dt4/runT4;
                    //Stage 5: Cross-boundary PMHL
                    stageDurations[PH2H_Cross]=intervalT-updateT;
                    dt5=stageDurations[PH2H_Cross]; qt5=dt5/runT5;
                }else{
                    dt4=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]; qt4=dt4/runT4;
                }
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }else{
            dt1=intervalT; qt1=dt1/runT1;
        }

        qtAll=qt1+qt2+qt4+qt5;
        cout<<"Stage 1 (Dijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (PCH): Duration: "<<dt2<<" ("<<stageDurations[PCH_No]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-Post): Duration: "<<dt4<<" ("<<stageDurations[PH2H_Post]<<") s; Throughput number: "<<qt4<<" ; query time: "<<1000 * runT4 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Extend): Duration: "<<dt5<<" ("<<stageDurations[PH2H_Cross]<<") s; Throughput number: "<<qt5<<" ; query time: "<<1000 * runT5 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
        stageQueryT[Dijk]+=runT1;
        stageQueryT[PCH_No]+=runT2;
        stageQueryT[PH2H_Post]+=runT4;
        stageQueryT[PH2H_Cross]+=runT5;
    }else if(algoChoice==6){
        runT1 = EffiVPLStage(ODpair,1000,Dijk,qTime1);
        runT2 = EffiVPLStage(ODpair, runtimes / 10, PCH_No, qTime2);
        runT4 = EffiVPLStage(ODpair, runtimes, PH2H_Post, qTime4);
        runT5 = EffiVPLStage(ODpair, runtimes, PH2H_Cross, qTime5);
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT) {
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT) {
                dt2=stageDurations[PCH_No]; qt2=dt2/runT2;
                //Stage 4: Post-boundary PMHL
                updateT+=stageDurations[PH2H_Post];
                if(updateT<intervalT) {
                    dt4=stageDurations[PH2H_Post]; qt4=dt4/runT4;
                    //Stage 5: Cross-boundary PMHL
                    stageDurations[PH2H_Cross]=intervalT-updateT;
                    dt5=stageDurations[PH2H_Cross]; qt5=dt5/runT5;
                }else{
                    dt4=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]; qt4=dt4/runT4;
                }
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }else{
            dt1=intervalT; qt1=dt1/runT1;
        }

        qtAll=qt1+qt2+qt4+qt5;
        cout<<"Stage 1 (Dijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (MPHL-CH): Duration: "<<dt2<<" ("<<stageDurations[PCH_No]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (MPHL-Post): Duration: "<<dt4<<" ("<<stageDurations[PH2H_Post]<<") s; Throughput number: "<<qt4<<" ; query time: "<<1000 * runT4 << " ms."<<endl;
        cout<<"Stage 4 (MPHL-Extend): Duration: "<<dt5<<" ("<<stageDurations[PH2H_Cross]<<") s; Throughput number: "<<qt5<<" ; query time: "<<1000 * runT5 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
        stageQueryT[Dijk]+=runT1;
        stageQueryT[PCH_No]+=runT2;
        stageQueryT[PH2H_Post]+=runT4;
        stageQueryT[PH2H_Cross]+=runT5;
    }else {
        cout<<"Wrong query type! "<<algoChoice<<endl; exit(1);
    }

}

vector<double> Graph::StageDurationCompute(int intervalT){
//    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    vector<double> durations(5,0);
    double updateT=0;//overall update time

    if(algoChoice==1){
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT){
            durations[Dijk]=stageDurations[Dijk];
            //Stage 2: CH
            updateT+=stageDurations[CH];
            if(updateT<intervalT){
                durations[CH]=stageDurations[CH];
                //Stage 3: CH
                stageDurations[H2H]=intervalT-updateT;
                durations[H2H]=stageDurations[H2H];
            }else{
                durations[CH]=intervalT-stageDurations[Dijk];
            }
        }else{
            durations[Dijk]=intervalT;
        }
//        stageUpdateT[Dijk]+=stageDurations[Dijk];
//        stageUpdateT[CH]+=stageDurations[CH];
//        stageUpdateT[H2H]+=stageDurations[H2H];
    }
    else if(algoChoice==3){
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT){
            durations[Dijk]=stageDurations[Dijk];
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT){
                durations[PCH_No]=stageDurations[PCH_No];
                //Stage 3: No-boundary PMHL
                updateT+=stageDurations[PH2H_No];
                if(updateT<intervalT){
                    durations[PH2H_No]=stageDurations[PH2H_No];
                    //Stage 4: Post-boundary PMHL
                    updateT+=stageDurations[PH2H_Post];
                    if(updateT<intervalT){
                        durations[PH2H_Post]=stageDurations[PH2H_Post];
                        //Stage 5: Cross-boundary PMHL
                        stageDurations[PH2H_Cross]=intervalT-updateT;
                        durations[PH2H_Cross]=stageDurations[PH2H_Cross];
                    }else{
                        durations[PH2H_Post]=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]-stageDurations[PH2H_No];
                    }
                }else{
                    durations[PH2H_No]=intervalT-stageDurations[Dijk]-stageDurations[PCH_No];
                }
            }else{
                durations[PCH_No]=intervalT-stageDurations[Dijk];
            }
        }
        else{
            durations[Dijk]=intervalT;
        }

//        stageUpdateT[Dijk]+=stageDurations[Dijk];
//        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
//        stageUpdateT[PH2H_No]+=stageDurations[PH2H_No];
//        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
//        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
    }
    else if(algoChoice==5){
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT) {
            durations[Dijk]=stageDurations[Dijk];
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT) {
                durations[PCH_No]=stageDurations[PCH_No];
                //Stage 4: Post-boundary PMHL
                updateT+=stageDurations[PH2H_Post];
                if(updateT<intervalT) {
                    durations[PH2H_Post]=stageDurations[PH2H_Post];
                    //Stage 5: Cross-boundary PMHL
                    stageDurations[PH2H_Cross]=intervalT-updateT;
                    durations[PH2H_Cross]=stageDurations[PH2H_Cross];
                }else{
                    durations[PH2H_Post]=intervalT-stageDurations[Dijk]-stageDurations[PCH_No];
                }
            }else{
                durations[PCH_No]=intervalT-stageDurations[Dijk];
            }
        }else{
            durations[Dijk]=intervalT;
        }

//        stageUpdateT[Dijk]+=stageDurations[Dijk];
//        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
//        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
//        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
    }
    else if(algoChoice==6){
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT) {
            durations[Dijk]=stageDurations[Dijk];
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT) {
                durations[PCH_No]=stageDurations[PCH_No];
                //Stage 4: Post-boundary PMHL
                updateT+=stageDurations[PH2H_Post];
                if(updateT<intervalT) {
                    durations[PH2H_Post]=stageDurations[PH2H_Post];
                    //Stage 5: Cross-boundary PMHL
                    stageDurations[PH2H_Cross]=intervalT-updateT;
                    durations[PH2H_Cross]=stageDurations[PH2H_Cross];
                }else{
                    durations[PH2H_Post]=intervalT-stageDurations[Dijk]-stageDurations[PCH_No];
                }
            }else{
                durations[PCH_No]=intervalT-stageDurations[Dijk];
            }
        }else{
            durations[Dijk]=intervalT;
        }

//        stageUpdateT[Dijk]+=stageDurations[Dijk];
//        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
//        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
//        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
    }else {
        cout<<"Wrong query type! "<<algoChoice<<endl; exit(1);
    }
    for(int i=0;i<durations.size();++i){
        if(i<=PCH_No){
            cout<<"Duration of Q-Stage "<<i+1<<" : "<<durations[i]<<" s"<<endl;
        }else{
            if((algoChoice==6||algoChoice==5) && i>=PH2H_Post){
                cout<<"Duration of Q-Stage "<<i<<" : "<<durations[i]<<" s"<<endl;
            }
        }
    }
    return durations;
}
//for VPL
void Graph::StageDurationComputeSchedule(int intervalT,vector<vector<double>>& durations){
//    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    durations.assign(VPLUpdateTimes.size(),vector<double>(5,0));

    for(int i=0;i<VPLUpdateTimes.size();++i){
        double updateT=0;//overall update time
        //Stage 1: Dijkstra
        updateT=VPLUpdateTimes[i][Dijk];
        if(updateT<intervalT) {
            durations[i][Dijk]=VPLUpdateTimes[i][Dijk];
            //Stage 2: PCH
            updateT+=VPLUpdateTimes[i][PCH_No];
            if(updateT<intervalT) {
                durations[i][PCH_No]=VPLUpdateTimes[i][PCH_No];
                //Stage 4: Post-boundary PMHL
                updateT+=VPLUpdateTimes[i][PH2H_Post];
                if(updateT<intervalT) {
                    durations[i][PH2H_Post]=VPLUpdateTimes[i][PH2H_Post];
                    //Stage 5: Cross-boundary PMHL
                    VPLUpdateTimes[i][PH2H_Cross]=intervalT-updateT;
                    durations[i][PH2H_Cross]=VPLUpdateTimes[i][PH2H_Cross];
                }else{
                    durations[i][PH2H_Post]=intervalT-VPLUpdateTimes[i][Dijk]-VPLUpdateTimes[i][PCH_No];
                }
            }else{
                durations[i][PCH_No]=intervalT-VPLUpdateTimes[i][Dijk];
            }
        }else{
            durations[i][Dijk]=intervalT;
        }
    }


//        stageUpdateT[Dijk]+=stageDurations[Dijk];
//        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
//        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
//        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];

}


//function for efficiency test
void Graph::EffiStageCheck(vector<pair<int,int>> & ODpair, int runtimes, vector<double> & queryTimes){
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;
    tt.start();
    double runT=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
    vector<double> qTime1, qTime2, qTime3, qTime4, qTime5;
    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    double updateT=0;//overall update time
    unsigned long long qt1=0,qt2=0,qt3=0,qt4=0,qt5=0,qtAll=0;
    int d1,d2;
    bool ifDebug=false;
//    ifDebug=true;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }
    clock_t start = clock();
    vector<int> results(runtimes,-1);

    if(algoChoice==1){
        runT1=EffiMHLStage(ODpair,runtimes/100,Dijk,qTime1);
        runT2=EffiMHLStage(ODpair,runtimes/10,CH,qTime2);
        runT3=EffiMHLStage(ODpair,runtimes,H2H,qTime3);
        queryTimes[Dijk]=runT1;
        queryTimes[CH]=runT2;
        queryTimes[H2H]=runT3;
    }
    else if(algoChoice==3){
        runT1=EffiPMHLStage(ODpair,runtimes/100,Dijk,qTime1);
        runT2=EffiPMHLStage(ODpair,runtimes/10,PCH_No,qTime2);
        runT3=EffiPMHLStage(ODpair,runtimes,PH2H_No,qTime3);
        runT4=EffiPMHLStage(ODpair,runtimes,PH2H_Post,qTime4);
        runT5=EffiPMHLStage(ODpair,runtimes,PH2H_Cross,qTime5);
        queryTimes[Dijk]=runT1;
        queryTimes[PCH_No]=runT2;
        queryTimes[PH2H_No]=runT3;
        queryTimes[PH2H_Post]=runT4;
        queryTimes[PH2H_Cross]=runT5;
    }
    else if(algoChoice==5){
        runT1 = EffiPostMHLStage(ODpair,runtimes/100,Dijk,qTime1);
        runT2 = EffiPostMHLStage(ODpair, runtimes/10, PCH_No,qTime2);
        runT4 = EffiPostMHLStage(ODpair, runtimes, PH2H_Post,qTime4);
        runT5 = EffiPostMHLStage(ODpair, runtimes, PH2H_Cross,qTime5);
        queryTimes[Dijk]=runT1;
        queryTimes[PCH_No]=runT2;
        queryTimes[PH2H_Post]=runT4;
        queryTimes[PH2H_Cross]=runT5;
    }else {
        cout<<"Wrong query type! "<<algoChoice<<endl; exit(1);
    }
    tt.stop();
    cout<<"Time for efficiency test: "<<tt.GetRuntime()<<" s."<<endl;
}

//function for efficiency test
void Graph::EffiStageCheck(vector<pair<int,int>> & ODpair, int runtimes, vector<vector<double>> & queryTimes){//return in seconds
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;
    tt.start();
    double runT=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
//    vector<double> qTime1, qTime2, qTime3, qTime4, qTime5;
    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    double updateT=0;//overall update time
    unsigned long long qt1=0,qt2=0,qt3=0,qt4=0,qt5=0,qtAll=0;
    int d1,d2;

    clock_t start = clock();
    vector<int> results(runtimes,-1);

    if(algoChoice==1){
        runT1=EffiMHLStage(ODpair,runtimes/100,Dijk,queryTimes[Dijk]);
        runT2=EffiMHLStage(ODpair,runtimes/10,CH,queryTimes[CH]);
        runT5=EffiMHLStage(ODpair,runtimes,H2H,queryTimes[H2H]);
//        aveQT[Dijk]=runT1;
//        aveQT[CH]=runT2;
//        aveQT[H2H]=runT5;
    }
    else if(algoChoice==3){
        runT1=EffiPMHLStage(ODpair,runtimes/100,Dijk,queryTimes[Dijk]);
        runT2=EffiPMHLStage(ODpair,runtimes/10,PCH_No,queryTimes[PCH_No]);
        runT3=EffiPMHLStage(ODpair,runtimes,PH2H_No,queryTimes[PH2H_No]);
        runT4=EffiPMHLStage(ODpair,runtimes,PH2H_Post,queryTimes[PH2H_Post]);
        runT5=EffiPMHLStage(ODpair,runtimes,PH2H_Cross,queryTimes[PH2H_Cross]);
//        aveQT[Dijk]=runT1;
//        aveQT[PCH_No]=runT2;
//        aveQT[PH2H_No]=runT3;
//        aveQT[PH2H_Post]=runT4;
//        aveQT[PH2H_Cross]=runT5;
    }
    else if(algoChoice==5){
        runT1 = EffiPostMHLStage(ODpair,runtimes/100,Dijk,queryTimes[Dijk]);
        runT2 = EffiPostMHLStage(ODpair, runtimes/10, PCH_No,queryTimes[PCH_No]);
        runT4 = EffiPostMHLStage(ODpair, runtimes, PH2H_Post,queryTimes[PH2H_Post]);
        runT5 = EffiPostMHLStage(ODpair, runtimes, PH2H_Cross,queryTimes[PH2H_Cross]);
//        aveQT[Dijk]=runT1;
//        aveQT[PCH_No]=runT2;
//        aveQT[PH2H_Post]=runT4;
//        aveQT[PH2H_Cross]=runT5;
    }
    else if(algoChoice==6){
        runT1 = EffiVPLStage(ODpair,runtimes/100,Dijk,queryTimes[Dijk]);
        runT2 = EffiVPLStage(ODpair, runtimes/10, PCH_No,queryTimes[PCH_No]);
        if(algoUpdate>=PH2H_Post){
            runT4 = EffiVPLStage(ODpair, runtimes, PH2H_Post,queryTimes[PH2H_Post]);
            if(algoUpdate==PH2H_Cross){
                runT5 = EffiVPLStage(ODpair, runtimes, PH2H_Cross,queryTimes[PH2H_Cross]);
            }
        }


//        aveQT[Dijk]=runT1;
//        aveQT[PCH_No]=runT2;
//        aveQT[PCH_Hybrid]=runT3;
//        aveQT[PH2H_Post]=runT4;
//        aveQT[PH2H_Cross]=runT5;
    }else {
        cout<<"Wrong query type! "<<algoChoice<<endl; exit(1);
    }
    tt.stop();
    cout<<"Time for efficiency test: "<<tt.GetRuntime()<<" s."<<endl;
//    return runT5;
}

//function for efficiency test
void Graph::EffiStageCheck(vector<pair<int,int>> & ODpair, int runtimes, vector<vector<double>> & queryTimes, vector<double> & aveQT){//return in seconds
    runtimes=min(runtimes,10000);
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;
    tt.start();
    double runT=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
//    vector<double> qTime1, qTime2, qTime3, qTime4, qTime5;
    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    double updateT=0;//overall update time
    unsigned long long qt1=0,qt2=0,qt3=0,qt4=0,qt5=0,qtAll=0;
    int d1,d2;

    clock_t start = clock();
    vector<int> results(runtimes,-1);

    if(algoChoice==1){
        runT1=EffiMHLStage(ODpair,min(runtimes,1000),Dijk,queryTimes[Dijk]);
        runT2=EffiMHLStage(ODpair,min(runtimes,1000),CH,queryTimes[CH]);
        runT5=EffiMHLStage(ODpair,runtimes,H2H,queryTimes[H2H]);
        aveQT[Dijk]=runT1;
        aveQT[CH]=runT2;
        aveQT[H2H]=runT5;
    }
    else if(algoChoice==3){
        runT1=EffiPMHLStage(ODpair,min(runtimes,1000),Dijk,queryTimes[Dijk]);
        runT2=EffiPMHLStage(ODpair,min(runtimes,1000),PCH_No,queryTimes[PCH_No]);
        runT3=EffiPMHLStage(ODpair,runtimes,PH2H_No,queryTimes[PH2H_No]);
        runT4=EffiPMHLStage(ODpair,runtimes,PH2H_Post,queryTimes[PH2H_Post]);
        runT5=EffiPMHLStage(ODpair,runtimes,PH2H_Cross,queryTimes[PH2H_Cross]);
        aveQT[Dijk]=runT1;
        aveQT[PCH_No]=runT2;
        aveQT[PH2H_No]=runT3;
        aveQT[PH2H_Post]=runT4;
        aveQT[PH2H_Cross]=runT5;
    }
    else if(algoChoice==5){
        runT1 = EffiPostMHLStage(ODpair, min(runtimes,1000),Dijk,queryTimes[Dijk]);
        runT2 = EffiPostMHLStage(ODpair, min(runtimes,1000), PCH_No,queryTimes[PCH_No]);
        runT4 = EffiPostMHLStage(ODpair, runtimes, PH2H_Post,queryTimes[PH2H_Post]);
        runT5 = EffiPostMHLStage(ODpair, runtimes, PH2H_Cross,queryTimes[PH2H_Cross]);
        aveQT[Dijk]=runT1;
        aveQT[PCH_No]=runT2;
        aveQT[PH2H_Post]=runT4;
        aveQT[PH2H_Cross]=runT5;
    }
    else if(algoChoice==6){
        runT1 = EffiVPLStage(ODpair,min(runtimes,1000),Dijk,queryTimes[Dijk]);
        runT2 = EffiVPLStage(ODpair, min(runtimes,1000), PCH_No,queryTimes[PCH_No]);
        if(algoUpdate>=PH2H_Post){
            runT4 = EffiVPLStage(ODpair, runtimes, PH2H_Post,queryTimes[PH2H_Post]);
            if(algoUpdate==PH2H_Cross){
                runT5 = EffiVPLStage(ODpair, runtimes, PH2H_Cross,queryTimes[PH2H_Cross]);
            }
        }

        aveQT[Dijk]=runT1;
        aveQT[PCH_No]=runT2;
        aveQT[PH2H_Post]=runT4;
        aveQT[PH2H_Cross]=runT5;
    }else {
        cout<<"Wrong query type! "<<algoChoice<<endl; exit(1);
    }
    tt.stop();
    cout<<"Time for efficiency test: "<<tt.GetRuntime()<<" s."<<endl;
//    return runT5;
}

//function for efficiency test
void Graph::EffiStageCheckReal(vector<pair<int,int>> & ODpair, int runtimes, vector<vector<double>> & queryTimes){//return in seconds

    int s, t;
    Timer tt;
    tt.start();
    double runT=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
//    vector<double> qTime1, qTime2, qTime3, qTime4, qTime5;
    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    double updateT=0;//overall update time
    unsigned long long qt1=0,qt2=0,qt3=0,qt4=0,qt5=0,qtAll=0;
    int d1,d2;

    clock_t start = clock();
    if(runtimes>ODpair.size()){
        runtimes=ODpair.size();
    }
    cout<<"Efficiency test. Run times: "<<runtimes<<" ("<<ODpair.size()<<")"<<endl;
    vector<int> results(runtimes,-1);

    if(algoChoice==1){
        runT1=EffiMHLStage(ODpair,min(runtimes,1000),Dijk,queryTimes[Dijk]);
        runT2=EffiMHLStage(ODpair,min(runtimes,1000),CH,queryTimes[CH]);
        runT5=EffiMHLStage(ODpair,runtimes,H2H,queryTimes[H2H]);
//        aveQT[Dijk]=runT1;
//        aveQT[CH]=runT2;
//        aveQT[H2H]=runT5;
    }
    else if(algoChoice==3){
        runT1=EffiPMHLStage(ODpair,min(runtimes,1000),Dijk,queryTimes[Dijk]);
        runT2=EffiPMHLStage(ODpair,min(runtimes,1000),PCH_No,queryTimes[PCH_No]);
        runT3=EffiPMHLStage(ODpair,runtimes,PH2H_No,queryTimes[PH2H_No]);
        runT4=EffiPMHLStage(ODpair,runtimes,PH2H_Post,queryTimes[PH2H_Post]);
        runT5=EffiPMHLStage(ODpair,runtimes,PH2H_Cross,queryTimes[PH2H_Cross]);
//        aveQT[Dijk]=runT1;
//        aveQT[PCH_No]=runT2;
//        aveQT[PH2H_No]=runT3;
//        aveQT[PH2H_Post]=runT4;
//        aveQT[PH2H_Cross]=runT5;
    }
    else if(algoChoice==5){
        runT1 = EffiPostMHLStage(ODpair, min(runtimes,1000),Dijk,queryTimes[Dijk]);
        runT2 = EffiPostMHLStage(ODpair, min(runtimes,1000), PCH_No,queryTimes[PCH_No]);
        runT4 = EffiPostMHLStage(ODpair, runtimes, PH2H_Post,queryTimes[PH2H_Post]);
        runT5 = EffiPostMHLStage(ODpair, runtimes, PH2H_Cross,queryTimes[PH2H_Cross]);
//        aveQT[Dijk]=runT1;
//        aveQT[PCH_No]=runT2;
//        aveQT[PH2H_Post]=runT4;
//        aveQT[PH2H_Cross]=runT5;
    }
    else if(algoChoice==6){
        runT1 = EffiVPLStage(ODpair, min(runtimes,1000),Dijk,queryTimes[Dijk]);
        runT2 = EffiVPLStage(ODpair, min(runtimes,1000), PCH_No,queryTimes[PCH_No]);
        if(algoUpdate>=PH2H_Post){
            runT4 = EffiVPLStage(ODpair, runtimes, PH2H_Post,queryTimes[PH2H_Post]);
            if(algoUpdate==PH2H_Cross){
                runT5 = EffiVPLStage(ODpair, runtimes, PH2H_Cross,queryTimes[PH2H_Cross]);
            }
        }

//        aveQT[Dijk]=runT1;
//        aveQT[PCH_No]=runT2;
//        aveQT[PCH_Hybrid]=runT3;
//        aveQT[PH2H_Post]=runT4;
//        aveQT[PH2H_Cross]=runT5;
    }else {
        cout<<"Wrong query type! "<<algoChoice<<endl; exit(1);
    }
    tt.stop();
    cout<<"Time for efficiency test: "<<tt.GetRuntime()<<" s."<<endl;
//    return runT5;
}

double Graph::analytical_update_first(double T_q, double T_u, double T_r, double T, double V_q) {//T is duration, T_u is the average update time for one update, return throughput
//	M/M/1 return min(1.0 / T_q - 1.0 / T_r, 1.0 / T_q - tau * T_u / (T * T_q));
//	M/G/1
    if (T_r < T_q)
        return 0;
    double a = 2 * (T_r - T_q) / (V_q + 2 * T_r * T_q - T_q * T_q);
    double lambda_u = T_u / T;
    double b = (1 - lambda_u) / T_q;
    double c = a < b ? a : b;
    if (c < 0) return 0;
    return c;
}

// function for estimating the system throughput, old
//double Graph::ThroughputEstimate(vector<vector<double>> &query_costs, vector<vector<double>> &update_costs, double threshold_time, double T) {
//    double queryT, duration, queryT_var;
//    double update_time;
//    double throughput=0.0;
//
//    for(int i=0;i<query_costs.size();++i){
//        duration = get_mean(update_costs[i]);
//        if(duration<=0)
//            continue;
//        queryT = get_mean(query_costs[i]);
//        queryT_var = get_var(query_costs[i]);
//        double tau = duration/T;
//        double thr = tau*analytical_update_first(queryT, 0, threshold_time/1000,  duration, queryT_var);
//        cout<<"Throughput of Q-stage "<<i+1<<" : "<<thr<<endl;
//        throughput+=thr;
//    }
//    return throughput;
//}

// function for estimating the system throughput
pair<double,double> Graph::ThroughputEstimate(vector<vector<double>> &query_costs, vector<vector<double>> &update_costs, double threshold_time, double T) {
    double queryT, duration, queryT_var;
    double update_time=0;
    double durationT=0;
    double throughput=0.0;
    double tau, thr;
    double query_time=INF;
//    vector<double> updateT;
//    updateT.push_back(0);

    for(int i=0;i<query_costs.size();++i){
        duration = get_mean(update_costs[i]);
//        updateT.push_back(update_time);
        if(duration<=0)
            continue;
        durationT+=duration;
        queryT = get_mean(query_costs[i]);
        queryT_var = get_var(query_costs[i]);
        tau = duration/T;
        if(queryT<query_time){
            query_time=queryT;
        }

        thr = analytical_update_first(queryT, update_time, threshold_time,  durationT, queryT_var);
//        cout<<i<<" . update time: "<<update_time<<" s ; estimate duration: "<<durationT<<" s"<<endl;
        update_time+=duration;
        cout<<"Throughput of Q-stage "<<i+1<<" : "<<thr<<" ; duration: "<< duration<<" s ; query time: "<<queryT*1000<<" ms"<<endl;
        throughput+=thr*tau;
    }
//    cout<<"Weighted average of throughput: "<<throughput<<endl;
    return make_pair(throughput,query_time);
}

//update first model
pair<double, double> Graph::simulator_UpdateFirst(vector<vector<double>> & query_cost, vector<vector<double>> & update_cost, vector<query> &queryList, int T, double period_time){
    int total_time; // terminal time for each period
    unsigned long long int count = 0; // the number of queries processed.

    double avg_response_time = 0.0;
    query* current_query;

    double c_time=0;  // current time for processing queries

//    T = T*microsecs_per_sec; // transform T to microseconds
    int num_updates = T/period_time;//number of batch updates
//    cout<<"batch number: "<<num_updates<<endl;
    int current_update_index=0;


    bool finished = false;

    int query_processed_in_one_update_slot = 0;

    double update_time_rest;


    while(current_update_index * period_time < T){

//        cout<<"batch "<<current_update_index<<endl;

        c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];
        double durationPoint=c_time;
        // start from c_time, end at total_time
        total_time = min((int)period_time * (current_update_index+1), T);
        query_processed_in_one_update_slot = 0;
        for(int i=0;i<query_cost.size();++i){
            update_time_rest = update_cost[i][current_update_index % update_cost.size()];//update time
//                cout<<"Duration of Q-Stage "<<i<<" : "<<update_time_rest<< " s"<<endl;
            c_time=durationPoint;
            if(update_time_rest==0)
                continue;

            durationPoint+=update_time_rest;
            while(c_time <= durationPoint) {//if there is remained time for querying
                if( count < queryList.size()) {
                    current_query = &queryList[count];
                }else {//if all queries are processed
                    finished = true;
//                    cout<<"!!! Finished the processing of all queries."<<endl;
                    break;
                }
//            cout<<"query "<<count<<": "<<current_query->init_time<<" "<<current_query->process_time<<" ";
                current_query->process_time = query_cost[i][query_processed_in_one_update_slot % query_cost[i].size()];//obtain the true query processing time
//            cout<<current_query->process_time<<endl;

                if (current_query->init_time> c_time){
                    c_time = current_query->init_time;
                }

                //c_time = max(c_time, current_query->init_time*1000000.0);
                if(c_time + current_query->process_time <= durationPoint) {
                    count++;
                    query_processed_in_one_update_slot++;
                    avg_response_time +=  (c_time -  current_query->init_time + current_query->process_time);//obtain the query response time
//                simulated_query_costs.push_back(current_query->process_time);
                    c_time += current_query->process_time;
                } else {//if the remained time is not enough for query processing
                    current_query->process_time = current_query->process_time - (durationPoint - c_time);
                    break;
                }
            }

        }
        if(finished) break;
//        if(c_time+1<total_time){
//            cout<<"Seems wrong. "<<c_time<<" "<<total_time<<" s"<<endl; exit(1);
//        }
        current_update_index++;
    }


    //cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
    return make_pair(count * 1.0 / T, avg_response_time / count);//return throughput (query per second) and average response time (in us)
}

pair<double,int> Graph::getFastestWorker(vector<double>& workers){
    pair<double,int> minT;
    minT.first=workers[0], minT.second=0;
    for(int i=1;i<workers.size();++i){
        if(workers[i]<minT.first){
            minT.first=workers[i]; minT.second=i;
        }
    }
    return minT;
}
//update first model, multiple workers
pair<double, double> Graph::simulator_UpdateFirst(vector<vector<double>> & query_cost, vector<vector<double>> & update_cost, vector<query> &queryList, int T, double period_time, int workerNum){
    int total_time; // terminal time for each period
    unsigned long long int count = 0; // the number of queries processed.

    double avg_response_time = 0.0;
    query* current_query;

    double c_time=0;  // current time for processing queries

//    T = T*microsecs_per_sec; // transform T to microseconds
    int num_updates = T/period_time;//number of batch updates
//    cout<<"batch number: "<<num_updates<<endl;
    int current_update_index=0;


    bool finished = false;

    int query_processed_in_one_update_slot = 0;

    double update_time_rest;

    if(workerNum>1){//with multiple workers
        benchmark::heap<2, int, long long int> pqueue(workerNum);
        int topID; long long int topValue;
        vector<double> workers(workerNum,0.0);
        vector<double> c_times(workerNum,0.0);
        pair<double, int> minT;
        while(current_update_index * period_time < T){
//        cout<<"batch "<<current_update_index<<endl;
            c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];
            double durationPoint=c_time;
            // start from c_time, end at total_time
            total_time = min((int)period_time * (current_update_index+1), T);
            query_processed_in_one_update_slot = 0;
            for(int i=0;i<query_cost.size();++i){
                update_time_rest = update_cost[i][current_update_index % update_cost.size()];//update time
//                cout<<"Duration of Q-Stage "<<i<<" : "<<update_time_rest<< " s"<<endl;
                c_time=durationPoint;
                for(int wi=0;wi<c_times.size();++wi){
                    c_times[wi]=c_time;
                    pqueue.update(wi,c_time*Resolution);
                }
                if(update_time_rest==0)
                    continue;

                durationPoint+=update_time_rest;
                while(c_time <= durationPoint && !pqueue.empty()) {//if there is remained time for querying

                    pqueue.extract_min(topID, topValue);
                    c_time=c_times[topID];
                    if( count < queryList.size()) {
                        current_query = &queryList[count];
                    }else {//if all queries are processed
                        finished = true;
//                    cout<<"!!! Finished the processing of all queries."<<endl;
                        break;
                    }

                    current_query->process_time = query_cost[i][query_processed_in_one_update_slot % query_cost[i].size()];//obtain the true query processing time
                    if (current_query->init_time> c_times[topID]){
                        c_times[topID] = current_query->init_time;
                    }

                    //c_time = max(c_time, current_query->init_time*1000000.0);
                    if(c_times[topID] + current_query->process_time <= durationPoint) {
                        count++;
                        query_processed_in_one_update_slot++;
                        avg_response_time +=  (c_times[topID] -  current_query->init_time + current_query->process_time);//obtain the query response time
//                simulated_query_costs.push_back(current_query->process_time);
                        c_times[topID] += current_query->process_time;
                        pqueue.update(topID,c_times[topID]*Resolution);
                    } else {//if the remained time is not enough for query processing
                        current_query->process_time = current_query->process_time - (durationPoint - c_time);
                        c_times[topID] += current_query->process_time;
                        break;
                    }
                    if(!pqueue.empty()){
                        c_time=c_times[pqueue.top_id()];
                    }else{
                        break;
                    }


                }

            }
            if(finished) break;
//        if(c_time+1<total_time){
//            cout<<"Seems wrong. "<<c_time<<" "<<total_time<<" s"<<endl; exit(1);
//        }
            current_update_index++;
        }
    }
    else{//single worker
        while(current_update_index * period_time < T){
//        cout<<"batch "<<current_update_index<<endl;
            c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];
            double durationPoint=c_time;
            // start from c_time, end at total_time
            total_time = min((int)period_time * (current_update_index+1), T);
            query_processed_in_one_update_slot = 0;
            for(int i=0;i<query_cost.size();++i){
                update_time_rest = update_cost[i][current_update_index % update_cost.size()];//update time
//                cout<<"Duration of Q-Stage "<<i<<" : "<<update_time_rest<< " s"<<endl;
                c_time=durationPoint;
                if(update_time_rest==0)
                    continue;

                durationPoint+=update_time_rest;
                while(c_time <= durationPoint) {//if there is remained time for querying
                    if( count < queryList.size()) {
                        current_query = &queryList[count];
                    }else {//if all queries are processed
                        finished = true;
//                    cout<<"!!! Finished the processing of all queries."<<endl;
                        break;
                    }
//            cout<<"query "<<count<<": "<<current_query->init_time<<" "<<current_query->process_time<<" ";
                    current_query->process_time = query_cost[i][query_processed_in_one_update_slot % query_cost[i].size()];//obtain the true query processing time
//            cout<<current_query->process_time<<endl;

                    if (current_query->init_time> c_time){
                        c_time = current_query->init_time;
                    }

                    //c_time = max(c_time, current_query->init_time*1000000.0);
                    if(c_time + current_query->process_time <= durationPoint) {
                        count++;
                        query_processed_in_one_update_slot++;
                        avg_response_time +=  (c_time -  current_query->init_time + current_query->process_time);//obtain the query response time
//                simulated_query_costs.push_back(current_query->process_time);
                        c_time += current_query->process_time;
                    } else {//if the remained time is not enough for query processing
                        current_query->process_time = current_query->process_time - (durationPoint - c_time);
                        break;
                    }
                }

            }
            if(finished) break;
//        if(c_time+1<total_time){
//            cout<<"Seems wrong. "<<c_time<<" "<<total_time<<" s"<<endl; exit(1);
//        }
            current_update_index++;
        }
    }




    //cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
    return make_pair(count * 1.0 / T, avg_response_time / count);//return throughput (query per second) and average response time (in us)
}
// true: affected, false: unaffected
bool Graph::judgeQueryAffect(int ID1, int ID2, vector<pair<bool,int>> & partiAffect, int Qstage_i){
    if(Qstage_i==PH2H_Cross) return false;
    int pid1=PartiTag[ID1].first, pid2=PartiTag[ID2].first;
    if(Qstage_i<=PCH_No){//is V-Stage 1
        if(pid1==pid2){
            if(!partiAffect[pid1].first){//unaffected
                innerBenefitNum++;
//                cout<<"V-Stage 1: unaffected same-partition query "<<ID1<<" "<<ID2<<endl;
                return false;
            }
        }
    }else if(Qstage_i<=PH2H_Post){//V-Stage 2
        if(pid1==pid2){
            if(partiAffect[pid1].second>AFF){
                innerBenefitNum++;
//                cout<<"V-Stage 2: unaffected same-partition query "<<ID1<<" "<<ID2<<endl;
                return false;
            }
        }else{//cross-partition query
            if(partiAffect[pid1].second==UNAFF && partiAffect[pid2].second==UNAFF){
                outerBenefitNum++;
//                cout<<"V-Stage 2: unaffected cross-partition query "<<ID1<<" "<<ID2<<endl;
                return false;
            }
        }
    }
    return true;
}

//update first model, multiple workers
pair<double, double> Graph::simulator_UpdateFirstVPL(vector<pair<int,int>>& ODpair, vector<vector<double>> & query_cost, vector<vector<double>> & update_cost, vector<query> &queryList, int T, double period_time, int workerNum){
    int total_time; // terminal time for each period
    unsigned long long int count = 0; // the number of queries processed.
    int ID1,ID2,query_i,query_j;
    double avg_response_time = 0.0;
    query* current_query;

    double c_time=0;  // current time for processing queries

//    T = T*microsecs_per_sec; // transform T to microseconds
    int num_updates = T/period_time;//number of batch updates
//    cout<<"batch number: "<<num_updates<<endl;
    int current_update_index=0;


    bool finished = false;

    int query_processed_in_one_update_slot = 0;

    double update_time_rest;

    if(workerNum>1){//with multiple workers
        benchmark::heap<2, int, long long int> pqueue(workerNum);
        int topID; long long int topValue;
        vector<double> workers(workerNum,0.0);
        vector<double> c_times(workerNum,0.0);
        pair<double, int> minT;
        while(current_update_index * period_time < T){
//        cout<<"batch "<<current_update_index<<endl;
            c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];
            double durationPoint=c_time;
            // start from c_time, end at total_time
            total_time = min((int)period_time * (current_update_index+1), T);
            query_processed_in_one_update_slot = 0;
            for(int i=0;i<query_cost.size();++i){
                update_time_rest = update_cost[i][current_update_index % update_cost.size()];//update time
//                cout<<"Duration of Q-Stage "<<i<<" : "<<update_time_rest<< " s"<<endl;
                c_time=durationPoint;
                for(int wi=0;wi<c_times.size();++wi){
                    c_times[wi]=c_time;
                    pqueue.update(wi,c_time*Resolution);
                }
                if(update_time_rest==0)
                    continue;

                durationPoint+=update_time_rest;
                while(c_time <= durationPoint && !pqueue.empty()) {//if there is remained time for querying

                    pqueue.extract_min(topID, topValue);
                    c_time=c_times[topID];
                    if( count < queryList.size()) {
                        current_query = &queryList[count];
                    }else {//if all queries are processed
                        finished = true;
//                    cout<<"!!! Finished the processing of all queries."<<endl;
                        break;
                    }

                    current_query->process_time = query_cost[i][query_processed_in_one_update_slot % query_cost[i].size()];//obtain the true query processing time
                    if (current_query->init_time> c_times[topID]){
                        c_times[topID] = current_query->init_time;
                    }

                    //c_time = max(c_time, current_query->init_time*1000000.0);
                    if(c_times[topID] + current_query->process_time <= durationPoint) {
                        count++;
                        query_processed_in_one_update_slot++;
                        avg_response_time +=  (c_times[topID] -  current_query->init_time + current_query->process_time);//obtain the query response time
//                simulated_query_costs.push_back(current_query->process_time);
                        c_times[topID] += current_query->process_time;
                        pqueue.update(topID,c_times[topID]*Resolution);
                    } else {//if the remained time is not enough for query processing
                        current_query->process_time = current_query->process_time - (durationPoint - c_time);
                        c_times[topID] += current_query->process_time;
                        break;
                    }
                    if(!pqueue.empty()){
                        c_time=c_times[pqueue.top_id()];
                    }else{
                        break;
                    }


                }

            }
            if(finished) break;
//        if(c_time+1<total_time){
//            cout<<"Seems wrong. "<<c_time<<" "<<total_time<<" s"<<endl; exit(1);
//        }
            current_update_index++;
        }
    }
    else{//single worker
        while(current_update_index * period_time < T){
//        cout<<"batch "<<current_update_index<<endl;
            c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];
            double durationPoint=c_time;
            // start from c_time, end at total_time
            total_time = min((int)period_time * (current_update_index+1), T);
            query_processed_in_one_update_slot = 0;
            for(int i=0;i<query_cost.size();++i){
                update_time_rest = update_cost[current_update_index][i];//update time
//                cout<<"Duration of Q-Stage "<<i<<" : "<<update_time_rest<< " s"<<endl;
                c_time=durationPoint;
                if(update_time_rest==0)
                    continue;

                durationPoint+=update_time_rest;
                while(c_time <= durationPoint) {//if there is remained time for querying
                    if( count < queryList.size()) {
                        current_query = &queryList[count];
                    }else {//if all queries are processed
                        finished = true;
//                    cout<<"!!! Finished the processing of all queries."<<endl;
                        break;
                    }
//            cout<<"query "<<count<<": "<<current_query->init_time<<" "<<current_query->process_time<<" ";
                    query_i=query_processed_in_one_update_slot % query_cost[i].size();
                    if(i<PH2H_Cross){
                        queryNum++;
                        query_j=query_processed_in_one_update_slot % ODpair.size();
                        ID1=ODpair[query_j].first, ID2=ODpair[query_j].second;
                        bool ifAff=judgeQueryAffect(ID1,ID2,partiAffectInfo[current_update_index],i);
                        if(ifAff){//if affected
                            current_query->process_time = query_cost[i][query_i];//obtain the true query processing time
                        }else{//if unaffected
                            current_query->process_time = query_cost[PH2H_Cross][query_i];//obtain the true query processing time
                        }
                    }else{
                        current_query->process_time = query_cost[i][query_i];//obtain the true query processing time
                    }


//            cout<<current_query->process_time<<endl;

                    if (current_query->init_time> c_time){
                        c_time = current_query->init_time;
                    }

                    //c_time = max(c_time, current_query->init_time*1000000.0);
                    if(c_time + current_query->process_time <= durationPoint) {
                        count++;
                        query_processed_in_one_update_slot++;
                        avg_response_time +=  (c_time -  current_query->init_time + current_query->process_time);//obtain the query response time
//                simulated_query_costs.push_back(current_query->process_time);
                        c_time += current_query->process_time;
                    } else {//if the remained time is not enough for query processing
                        current_query->process_time = current_query->process_time - (durationPoint - c_time);
                        break;
                    }
                }

            }
            if(finished) break;
//        if(c_time+1<total_time){
//            cout<<"Seems wrong. "<<c_time<<" "<<total_time<<" s"<<endl; exit(1);
//        }
            current_update_index++;
        }
    }



//    queryNum=count;
    //cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
    return make_pair(count * 1.0 / T, avg_response_time / count);//return throughput (query per second) and average response time (in us)
}


// function for simulating the system throughput
double Graph::ThroughputSimulate(vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime, int workerNum){//period_time is update interval time, simulation_time is the overall time for simulation
    Timer tt;
    tt.start();
    double l = 0, r = 2000000; // be careful of 'r'
    double throughput = 0.0;
    int lambda;
    int gap=1;

    r*=workerNum;

    if(r>workerNum*10/queryTime){
        r=workerNum*10/queryTime;
    }

    cout<<"Simulate the throughput of system... r: "<<r<<" ; gap: "<<gap<<" ; worker number: "<< workerNum<<endl;
    while(l < r){//why use different lambda?
        lambda = (l+r)/2;
        if(lambda==0){
            lambda=1;
        }
        vector<query> queryList;
        lambdaCount[lambda]++;
        if(lambdaCache.find(lambda)==lambdaCache.end()){//if not found
            queryList = generate_queryList(lambda, simulation_time);
            if(lambdaCache.size()<=10){//if cache is enough
                lambdaCache.insert({lambda,queryList});
                lambdaPQueue.update(lambda,lambdaCount[lambda]);
            }else{//if cache is not enough
                if(lambdaCount[lambda]>lambdaCount[lambdaPQueue.top_id()]){//if this is more frequent than the top one in pqueue
                    int topID, topCount;
                    lambdaPQueue.extract_min(topID,topCount);
                    cout<<"Replace "<<topID<<"("<<topCount<<") with "<<lambda<<"("<<lambdaCount[lambda]<<"). ";
                    lambdaCache.erase(topID);
                    lambdaCache.insert({lambda,queryList});
                    lambdaPQueue.update(lambda,lambdaCount[lambda]);
                }
            }
        }else{//if found in cache
            cout<<"Find in cache. ";
            lambdaPQueue.update(lambda,lambdaCount[lambda]);
            queryList=lambdaCache[lambda];
        }
//        auto queryList = generate_queryList(lambda, simulation_time);
        cout<<"lambda: "<<lambda<<" ; the size of querylist: "<<queryList.size();

//        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time);
        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time, workerNum);

//        pair<double,double> result = simulator_UpdateFirst(simulated_query_costs, simulated_update_costs, queryList, simulation_time, period_time);
        if(result.second <= threshold_time){
            if(throughput < result.first) {
//                cout<<"here:" << result.first<<endl;
                throughput=result.first;

                cout<<" ; throughput: "<<result.first<< " ; response query time: "<< result.second*1000 <<" ms";
            }
            l = lambda + gap;
//            l=l*1.1;
        }else{
            r = lambda - gap;
//            r=r/1.1;
        }
        cout<<endl;
    }

    tt.stop();
    cout<<"Time for throughput simulation: "<<tt.GetRuntime()<<" s"<<endl;
    return throughput;

}

double Graph::ThroughputSimulateReal(vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime, int workerNum){//period_time is update interval time, simulation_time is the overall time for simulation
    Timer tt;
    tt.start();
    double l = 0, r = 2000000; // be careful of 'r'
    double throughput = 0.0;
    int lambda;
    int gap=1;

    r*=workerNum;

    if(r>workerNum*4/queryTime){
        r=workerNum*4/queryTime;
    }

    cout<<"Simulate the throughput of system... r: "<<r<<" ; gap: "<<gap<<" ; worker number: "<< workerNum<<endl;
    while(l < r){//why use different lambda?
        lambda = (l+r)/2;
        if(lambda==0){
            lambda=1;
        }
        vector<query> queryList;
        lambdaCount[lambda]++;
        if(lambdaCache.find(lambda)==lambdaCache.end()){//if not found
            queryList = generate_queryList(lambda, simulation_time);
            if(lambdaCache.size()<=lambdaCacheSize){//if cache is enough
                lambdaCache.insert({lambda,queryList});
                lambdaPQueue.update(lambda,lambdaCount[lambda]);
            }else{//if cache is not enough
                if(lambdaCount[lambda]>lambdaCount[lambdaPQueue.top_id()]){//if this is more frequent than the top one in pqueue
                    int topID, topCount;
                    lambdaPQueue.extract_min(topID,topCount);
                    cout<<"Replace "<<topID<<"("<<topCount<<") with "<<lambda<<"("<<lambdaCount[lambda]<<"). ";
                    lambdaCache.erase(topID);
                    lambdaCache.insert({lambda,queryList});
                    lambdaPQueue.update(lambda,lambdaCount[lambda]);
                }
            }
        }else{//if found in cache
            cout<<"Find in cache. ";
            lambdaPQueue.update(lambda,lambdaCount[lambda]);
            queryList=lambdaCache[lambda];
        }
//        auto queryList = generate_queryList(lambda, simulation_time);
        cout<<"lambda: "<<lambda<<" ; the size of querylist: "<<queryList.size();

//        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time);
        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time, workerNum);

//        pair<double,double> result = simulator_UpdateFirst(simulated_query_costs, simulated_update_costs, queryList, simulation_time, period_time);
        if(result.second <= threshold_time){
            if(throughput < result.first) {
//                cout<<"here:" << result.first<<endl;
                throughput=result.first;
//                cout<<"lambda: "<<lambda<<" ; the size of querylist: "<<queryList.size();
                cout<<" ; throughput: "<<result.first<< " ; response query time: "<< result.second*1000 <<" ms";
            }
            l = lambda + gap;
//            l=l*1.1;
        }else{
            r = lambda - gap;
//            r=r/1.1;
        }
        cout<<endl;
    }

    tt.stop();
    cout<<"Throughput: "<<throughput<<" . Time for throughput simulation: "<<tt.GetRuntime()<<" s"<<endl;
    return throughput;

}

// function for simulating the system throughput
double Graph::ThroughputSimulateRealVPL(vector<pair<int,int>>& ODpair,vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime, int workerNum){//period_time is update interval time, simulation_time is the overall time for simulation
    Timer tt;
    tt.start();
    double l = 0, r = 2000000; // be careful of 'r'
    double throughput = 0.0;
    int lambda;
    int gap=1;

    r*=workerNum;

    if(r>workerNum*4/queryTime){
        r=workerNum*4/queryTime;
    }

    cout<<"Simulate the throughput of system... r: "<<r<<" ; gap: "<<gap<<" ; worker number: "<< workerNum<<endl;
    while(l < r){//why use different lambda?
        lambda = (l+r)/2;
        if(lambda<=0){
            lambda=1;
        }
        vector<query> queryList;
        lambdaCount[lambda]++;
        if(lambdaCache.find(lambda)==lambdaCache.end()){//if not found
            queryList = generate_queryList(lambda, simulation_time);
            if(lambdaCache.size()<=lambdaCacheSize){//if cache is enough
                lambdaCache.insert({lambda,queryList});
                lambdaPQueue.update(lambda,lambdaCount[lambda]);
            }else{//if cache is not enough
                if(lambdaCount[lambda]>lambdaCount[lambdaPQueue.top_id()]){//if this is more frequent than the top one in pqueue
                    int topID, topCount;
                    lambdaPQueue.extract_min(topID,topCount);
                    cout<<"Replace "<<topID<<"("<<topCount<<") with "<<lambda<<"("<<lambdaCount[lambda]<<"). ";
                    lambdaCache.erase(topID);
                    lambdaCache.insert({lambda,queryList});
                    lambdaPQueue.update(lambda,lambdaCount[lambda]);
                }
            }
        }else{//if found in cache
            cout<<"Find in cache. ";
            lambdaPQueue.update(lambda,lambdaCount[lambda]);
            queryList=lambdaCache[lambda];
        }

        cout<<"lambda: "<<lambda<<" ; the size of querylist: "<<queryList.size();
        innerBenefitNum=0, outerBenefitNum=0; queryNum=0;
//        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time);
        pair<double,double> result = simulator_UpdateFirstVPL(ODpair, query_costs, update_costs, queryList, simulation_time, period_time, workerNum);

//        pair<double,double> result = simulator_UpdateFirst(simulated_query_costs, simulated_update_costs, queryList, simulation_time, period_time);
        if(result.second <= threshold_time){
            if(throughput < result.first) {
//                cout<<"here:" << result.first<<endl;
                throughput=result.first;
//                cout<<"lambda: "<<lambda<<" ; the size of querylist: "<<queryList.size();
                cout<<" ; throughput: "<<result.first<< " ; response query time: "<< result.second*1000 <<" ms ; unaffected query number: "<<innerBenefitNum+outerBenefitNum<<" ( "<<queryNum<<" )";
            }
//            gap=(r-l)/10;
//            if(gap<5) gap=1;
            l = lambda + gap;
//            l=l*1.1;
        }else{
//            gap=(r-l)/10;
//            if(gap<5) gap=1;
            r = lambda - gap;
//            r=r/1.1;
        }
        cout<<endl;
    }
    tt.stop();
    unaffectedQNum.push_back(innerBenefitNum+outerBenefitNum);
    totalQNum.push_back(queryNum);

    cout<<"Throughput: "<<throughput<<" . Time for throughput simulation: "<<tt.GetRuntime()<<" s"<<endl;
    return throughput;

}

//function for efficiency test
void Graph::GetBatchThroughput(vector<double> & queryTimes, int intervalT, unsigned long long & throughputNum, vector<double>& stageUpdateT){

    int s, t;
    Timer tt;

    double runT=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    double updateT=0;//overall update time
    unsigned long long qt1=0,qt2=0,qt3=0,qt4=0,qt5=0,qtAll=0;
    int d1,d2;
    bool ifDebug=false;
//    ifDebug=true;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }
    clock_t start = clock();

    if(algoChoice==1){
        runT1=queryTimes[Dijk];
        runT2=queryTimes[CH];
        runT3=queryTimes[H2H];
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT){
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: CH
            updateT+=stageDurations[CH];
            if(updateT<intervalT){
                dt2=stageDurations[CH]; qt2=dt2/runT2;
                //Stage 3: CH
                stageDurations[H2H]=intervalT-updateT;
                dt3=stageDurations[H2H]; qt3=dt3/runT3;
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }else{
            dt1=intervalT; qt1=dt1/runT1;
        }


        qtAll=qt1+qt2+qt3;
        cout<<"Stage 1 (BiDijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (CH): Duration: "<<dt2<<" ("<<stageDurations[CH]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (H2H): Duration: "<<dt3<<" ("<<stageDurations[H2H]<<") s; Throughput number: "<<qt3<<" ; query time: "<<1000 * runT3 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[CH]+=stageDurations[CH];
        stageUpdateT[H2H]+=stageDurations[H2H];
    }
    else if(algoChoice==3){
        runT1=queryTimes[Dijk];
        runT2=queryTimes[PCH_No];
        runT3=queryTimes[PH2H_No];
        runT4=queryTimes[PH2H_Post];
        runT5=queryTimes[PH2H_Cross];
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT){
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT){
                dt2=stageDurations[PCH_No]; qt2=dt2/runT2;
                //Stage 3: No-boundary PMHL
                updateT+=stageDurations[PH2H_No];
                if(updateT<intervalT){
                    dt3=stageDurations[PH2H_No]; qt3=dt3/runT3;
                    //Stage 4: Post-boundary PMHL
                    updateT+=stageDurations[PH2H_Post];
                    if(updateT<intervalT){
                        dt4=stageDurations[PH2H_Post]; qt4=dt4/runT4;
                        //Stage 5: Cross-boundary PMHL
                        stageDurations[PH2H_Cross]=intervalT-updateT;
                        dt5=stageDurations[PH2H_Cross]; qt5=dt5/runT5;
                    }else{
                        dt4=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]-stageDurations[PH2H_No]; qt4=dt4/runT4;
                    }
                }else{
                    dt3=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]; qt3=dt3/runT3;
                }
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }
        else{
            dt1=intervalT; qt1=dt1/runT1;
        }

        qtAll=qt1+qt2+qt3+qt4+qt5;
        cout<<"Stage 1 (Dijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (PCH): Duration: "<<dt2<<" ("<<stageDurations[PCH_No]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-No): Duration: "<<dt3<<" ("<<stageDurations[PH2H_No]<<") s; Throughput number: "<<qt3<<" ; query time: "<<1000 * runT3 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Post): Duration: "<<dt4<<" ("<<stageDurations[PH2H_Post]<<") s; Throughput number: "<<qt4<<" ; query time: "<<1000 * runT4 << " ms."<<endl;
        cout<<"Stage 5 (PH2H-Extend): Duration: "<<dt5<<" ("<<stageDurations[PH2H_Cross]<<") s; Throughput number: "<<qt5<<" ; query time: "<<1000 * runT5 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
        stageUpdateT[PH2H_No]+=stageDurations[PH2H_No];
        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
    }
    else if(algoChoice==5){
        runT1 = queryTimes[Dijk];
        runT2 = queryTimes[PCH_No];
        runT4 = queryTimes[PH2H_Post];
        runT5 = queryTimes[PH2H_Cross];
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT) {
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT) {
                dt2=stageDurations[PCH_No]; qt2=dt2/runT2;
                //Stage 4: Post-boundary PMHL
                updateT+=stageDurations[PH2H_Post];
                if(updateT<intervalT) {
                    dt4=stageDurations[PH2H_Post]; qt4=dt4/runT4;
                    //Stage 5: Cross-boundary PMHL
                    stageDurations[PH2H_Cross]=intervalT-updateT;
                    dt5=stageDurations[PH2H_Cross]; qt5=dt5/runT5;
                }else{
                    dt4=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]; qt4=dt4/runT4;
                }
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }else{
            dt1=intervalT; qt1=dt1/runT1;
        }

        qtAll=qt1+qt2+qt4+qt5;
        cout<<"Stage 1 (Dijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (PCH): Duration: "<<dt2<<" ("<<stageDurations[PCH_No]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-Post): Duration: "<<dt4<<" ("<<stageDurations[PH2H_Post]<<") s; Throughput number: "<<qt4<<" ; query time: "<<1000 * runT4 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Extend): Duration: "<<dt5<<" ("<<stageDurations[PH2H_Cross]<<") s; Throughput number: "<<qt5<<" ; query time: "<<1000 * runT5 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
    }else {
        cout<<"Wrong query type! "<<algoChoice<<endl; exit(1);
    }

}

double Graph::EffiMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType, vector<double> &qTime) {
    algoQuery=queryType;
    Timer tt;
    int ID1, ID2, d1, d2;
    double runT=0;
    vector<int> results(runtimes,0);
    bool ifDebug=false;
//    ifDebug=true;

    for(int i=0;i<runtimes;++i){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
        tt.start();
        d1=QueryNP(ID1,ID2);
        results[i]=d1;
        tt.stop();
        runT+=tt.GetRuntime();
        qTime.push_back(tt.GetRuntime());
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                cout<<"Wrong result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
            }
        }
    }


    cout<<"Average Efficiency of Stage "<<queryType<<" : "<<1000*runT/runtimes<<" ms."<<endl;
    return runT/runtimes;
}

double Graph::EffiPMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType, vector<double> &qTime) {
    algoQuery=queryType;
    Timer tt;
    int ID1, ID2, d1, d2;
    double runT=0;
    vector<int> results(runtimes,0);
    bool ifDebug=false;
//    ifDebug=true;

    for(int i=0;i<runtimes;++i){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
        tt.start();
        d1=QueryPMHL(ID1,ID2);
        results[i]=d1;
        tt.stop();
        runT+=tt.GetRuntime();
        qTime.push_back(tt.GetRuntime());
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                cout<<"Wrong result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
            }
        }
    }

    cout<<"Average Efficiency of Stage "<<queryType<<" : "<<1000*runT/runtimes<<" ms."<<endl;
    return runT/runtimes;
}

double Graph::EffiPostMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType, vector<double> &qTime) {
    algoQuery=queryType;
    Timer tt;
    int ID1, ID2, d1, d2;
    double runT=0;
    vector<int> results(runtimes,0);
    bool ifDebug=false;
//    ifDebug=true;

    for(int i=0;i<runtimes;++i){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
        tt.start();
        d1=QueryPostMHL(ID1,ID2);
        results[i]=d1;
        tt.stop();
        runT+=tt.GetRuntime();
        qTime.push_back(tt.GetRuntime());
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                cout<<"Wrong result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
            }
        }
    }

    cout<<"Average Efficiency of Stage "<<queryType<<" : "<<1000*runT/runtimes<<" ms."<<endl;
    return runT/runtimes;
}

double Graph::EffiVPLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType, vector<double> &qTime) {
    algoQuery=queryType;
    Timer tt;
    int ID1, ID2, d1, d2;
    double runT=0;
    vector<int> results(runtimes,0);
    bool ifDebug=false;
//    ifDebug=true;

    for(int i=0;i<runtimes;++i){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
        tt.start();
        d1=QueryVPL(ID1,ID2);
        results[i]=d1;
        tt.stop();
        runT+=tt.GetRuntime();
        qTime.push_back(tt.GetRuntime());
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                cout<<"Wrong result! "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<") "<<d1<<" "<<d2<<endl; exit(1);
            }
        }
    }

    cout<<"Average Efficiency of Stage "<<queryType+1<<" : "<<1000*runT/runtimes<<" ms."<<endl;
    return runT/runtimes;
}

void Graph::DFSTree(vector<int>& tNodes, int id){
    tNodes.push_back(Tree[id].uniqueVertex);
    for(int i=0;i<Tree[id].ch.size();++i){
        DFSTree(tNodes,Tree[id].ch[i]);
    }
}





/// Index Maintenance
//function of testing the throughput on real-life updates
void Graph::RealUpdateThroughputTest(string updateFile){
    bool ifDebug=false;
    ifDebug=true;
    int runtimes=10000;
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    string line;
    vector<string> vs;
    int ID1,ID2,oldW,newW,weight;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==2);
    int batchNum=stoi(vs[0]);
    int batchInterval=stoi(vs[1]);
    int batchSize;
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNum);
    long long int aveBatchSize=0;
    int maxBatchSize=0, minBatchSize=INT32_MAX;
    for(int i=0;i<batchNum;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "));
        batchSize=stoi(vs[0]);
        assert(vs.size()==3*batchSize+1);
        aveBatchSize+=batchSize;
        if(maxBatchSize<batchSize) maxBatchSize=batchSize;
        if(minBatchSize>batchSize) minBatchSize=batchSize;
        for(int j=0;j<batchSize;++j){
            ID1=stoi(vs[3*j+1]), ID2=stoi(vs[3*j+2]), weight=stoi(vs[3*j+3]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    IF.close();
    aveBatchSize/=batchNum;
    cout<<"Update batch number: "<<batchNum<<" ; Average batch size: "<<aveBatchSize<<" ; Maximal batch size: "<<maxBatchSize<<" ; Minimal batch size: "<<minBatchSize<<" ; Batch interval: "<< batchInterval<<endl;

    string queryF = sourcePath+"/"+dataset + ".query";
    if(samePartiPortion!=-1){
        queryF=sourcePath+"/tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".sameParti_"+to_string(samePartiPortion)+".query";
    }
    ifstream IF2(queryF);
    if(!IF2){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF2>>num;
    for(int k=0;k<num;k++){
        IF2>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF2.close();

    //index maintenance
    Timer tt;
    Timer tRecord;
    double runT1=0, runT2 = 0;
    unsigned long long throughputNum=0;
    if(algoChoice==5){
        ProBeginVertexSetParti.assign(partiNum,vector<int>());
        vertexIDChLParti.assign(partiNum,set<int>());
        ProBeginVertexSetPartiExtend.assign(partiNum,vector<int>());
    }
    double updateTime=0;
    vector<double> stageUpdateT(5,0);
    vector<double> stageQueryT(5,0);
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;
    stageDurations.assign(5,0);
    EffiStageCheck(ODpair, runtimes, stageQueryT);
//    batchNum=2;
    for(int i=0;i<batchNum;++i){
        wBatchDec.clear(); wBatchInc.clear();
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[i].size();++j){
            ID1=batchUpdates[i][j].first.first, ID2=batchUpdates[i][j].first.second, weight=batchUpdates[i][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                if(it->first==ID2){
                    ifFind=true;
                    oldW=it->second;
                    if(oldW>weight){
//                        cout<<"Dec "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                        wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                    }else if(it->second<weight){
//                        cout<<"Inc "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                        wBatchInc.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                    }
                    break;
                }
            }
            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),weight});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<weight<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<i<<" . Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;

        int id1=2, id2=3;
//        cout<<"Before update "<<i<<". "<<id1<<"("<<NodeOrder[id1]<<","<<PartiTags[id1].first<<") "<<id2<<"("<<NodeOrder[id2]<<","<<PartiTags[id2].first<<") "<<QueryPostMHL(id1,id2)<<"("<<Dijkstra(id1,id2,Neighbor)<<")"<<endl;
//        cout<<"Before update "<<i<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") "<<endl;
//        for(auto it=Tree[rank[id1]].vert.begin();it!=Tree[rank[id1]].vert.end();++it){
//            if(it->first==id2){
//                cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//                break;
//            }
//        }

        //Step 1: Decrease updates
        tt.start();
        if(!wBatchDec.empty()){
            cout<<"Decrease update. "<<wBatchDec.size()<<endl;
            tRecord.start();
//                    boost::thread_group thread;

            if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::DecBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT1)));
                DecBatchThroughputNP(wBatchDec, i, runT1);
            }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                PMHLBatchUpdateDec(wBatchDec, i, runT1);
            }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                PostMHLBatchUpdateDec(wBatchDec, i, runT1);
            }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

//            if(ifDebug){
//                CorrectnessCheck(100);
//            }
        }

        //Step 2: Increase updates
        if(!wBatchInc.empty()){
            cout<<"Increase update. "<<wBatchInc.size()<<endl;
//            tRecord.start();
////            boost::thread_group thread;
            if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::IncBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT2)));
                IncBatchThroughputNP(wBatchInc, i, runT2);
            }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));
                PMHLBatchUpdateInc(wBatchInc, i, runT2);
            }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));//increase update is incorrect
                PostMHLBatchUpdateInc(wBatchInc, i, runT2);
            }

//            if(i==19){
//                ofstream OF(sourcePath+"/"+dataset+".CH19");
//                if(!OF.is_open()){
//                    cout<<"Fail to open file"<<endl; exit(1);
//                }
//                OF<<node_num<<" "<<edge_num<<endl;
//                for(int i=0;i<Neighbor.size();++i){
//                    for(auto it=Neighbor[i].begin();it!=Neighbor[i].end();++it){
//                        OF<<i<<" "<<it->first<<" "<<it->second<<endl;
//                    }
//                }
//                OF.close();
//                PostMHLIndexStoreCH(sourcePath+"/"+dataset+".CHIndex19");
//                cout<<"Write done."<<endl;
//            }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

//                    cout<<"After update: "<<endl;
//                    for(int i=0;i<Tree[rank[t]].disPost.size();++i){
//                        if(Tree[rank[t]].vAncestor[i]==hub){
//                            cout<<"Cnt: "<<Tree[rank[t]].disPost[i]<<" "<<Tree[rank[t]].cntPost[i]<<endl;
//                        }
//                    }
//                    for(auto it=Tree[rank[t]].vert.begin();it!=Tree[rank[t]].vert.end();++it){
//                        cout<<"sc: "<<t<<" "<<it->first<<" "<<it->second.first<<"("<<it->second.second<<")"<<endl;
//                    }
//                    cout<<"disB: "<<s<<" "<<bv<<" "<<Tree[rank[s]].disInf[BoundVertexMap[3][bv]]<<endl;


        }

        tt.stop();
        updateTime+=tt.GetRuntime();
        cout<<"Batch "<<i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;

        if(ifDebug){
            CorrectnessCheck(100);
        }

//        EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
        GetBatchThroughput(stageQueryT,batchInterval,throughputNum,stageUpdateT);
    }
//    for(int i=0;i<stageDurations.size();++i){
//        stageDurations[i]/=batchNum;
//    }
//    EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
    AverageStagePerformance(batchNum,stageUpdateT,stageQueryT);
    cout<<"\nPartiNum: "<<partiNum<<". Overall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average batch update Time: "<<updateTime/batchNum<<" s."<<endl;
}

void Graph::PostMHLIndexStoreCH(string filename){
    ofstream OF(filename);
    if(!OF.is_open()){
        cout<<"Cannot open file. "<<filename<<endl; exit(1);
    }
    cout<<"Writing labels..."<<endl;
    OF<<Tree.size()<<endl;
    for(int i=0;i<Tree.size();++i){
        OF<<Tree[i].uniqueVertex<<" "<<Tree[i].vert.size();
        for(int j=0;j!=Tree[i].vert.size();++j){
            OF<<" "<<Tree[i].vert[j].first<<" "<<Tree[i].vert[j].second.first<<" "<<Tree[i].vert[j].second.second;
        }
        OF<<endl;
    }
    OF.close();
}

void Graph::MHLIndexCompareCH(string filename){
    vector<map<int,pair<int,int>>> CHIndex;
    ifstream IF(filename);
    if(!IF.is_open()){
        cout<<"Cannot open file. "<<filename<<endl; exit(1);
    }
//    ofstream OF(filename+".result");
//    if(!OF.is_open()){
//        cout<<"Cannot open file. "<<filename+".result"<<endl; exit(1);
//    }
    cout<<"Index file: "<<filename<<endl;
    cout<<"Read and compare..."<<endl;
    string line;
    vector<string> vs;
    getline(IF,line);
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    int vNum=stoi(vs[0]);
    if(vNum!=node_num) {
        cout<<"Wrong node number. "<<vNum<<" "<<node_num<<endl; exit(1);
    }
    CHIndex.assign(vNum,map<int,pair<int,int>>());
    int ID1,ID2,dis,num,cnt;
    for(int i=0;i<vNum;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]), num=stoi(vs[1]);
        assert(vs.size()==3*num+2);
        for(int j=0;j<num;++j){
            ID2=stoi(vs[3*j+2]), dis=stoi(vs[3*j+3]), cnt=stoi(vs[3*j+4]);
            CHIndex[ID1].insert({ID2,make_pair(dis,cnt)});
        }
    }
    IF.close();
    for(int i=0;i<Tree.size();++i){
        ID1=Tree[i].uniqueVertex;
        for(auto it=Tree[i].vert.begin();it!=Tree[i].vert.end();++it){
            ID2=it->first, dis=it->second.first, cnt=it->second.second;
            if(CHIndex[ID1].find(ID2)!=CHIndex[ID1].end()){//if found
                if(CHIndex[ID1][ID2].first!=dis){
                    cout<<"Wrong label! "<<ID1<<"("<<NodeOrder[ID1]<<") "<<ID2<<"("<<NodeOrder[ID2]<<") "<<CHIndex[ID1][ID2].first<<"("<<dis<<") "<<CHIndex[ID1][ID2].second<<"("<<cnt<<")"<<endl;
//                    exit(1);
                }else if(CHIndex[ID1][ID2].second!=cnt ){
                    cout<<"Wrong cnt. "<<ID1<<"("<<NodeOrder[ID1]<<") "<<ID2<<"("<<NodeOrder[ID2]<<") "<<CHIndex[ID1][ID2].first<<"("<<dis<<") "<<CHIndex[ID1][ID2].second<<"("<<cnt<<")"<<endl;
                }
            }else{
                cout<<"Wrong. Not found. "<<ID1<<" "<<ID2<<" "<<dis<<endl; exit(1);
            }
        }
    }
}

void Graph::PostMHLIndexCompareCH(string filename){
    vector<map<int,pair<int,int>>> CHIndex;
    ifstream IF(filename);
    if(!IF.is_open()){
        cout<<"Cannot open file. "<<filename<<endl; exit(1);
    }
//    ofstream OF(filename+".result");
//    if(!OF.is_open()){
//        cout<<"Cannot open file. "<<filename+".result"<<endl; exit(1);
//    }
    cout<<"Index file: "<<filename<<endl;
    cout<<"Read and compare..."<<endl;
    string line;
    vector<string> vs;
    getline(IF,line);
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    int vNum=stoi(vs[0]);
    if(vNum!=node_num) {
        cout<<"Wrong node number. "<<vNum<<" "<<node_num<<endl; exit(1);
    }
    CHIndex.assign(vNum,map<int,pair<int,int>>());
    int ID1,ID2,dis,num,cnt;
    for(int i=0;i<vNum;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]), num=stoi(vs[1]);
        assert(vs.size()==3*num+2);
        for(int j=0;j<num;++j){
            ID2=stoi(vs[3*j+2]), dis=stoi(vs[3*j+3]), cnt=stoi(vs[3*j+4]);
            CHIndex[ID1].insert({ID2,make_pair(dis,cnt)});
        }
    }
    IF.close();
    for(int i=0;i<Tree.size();++i){
        ID1=Tree[i].uniqueVertex;
        for(auto it=Tree[i].vert.begin();it!=Tree[i].vert.end();++it){
            ID2=it->first, dis=it->second.first, cnt=it->second.second;
            if(CHIndex[ID1].find(ID2)!=CHIndex[ID1].end()){//if found
                if(CHIndex[ID1][ID2].first!=dis){
                    cout<<"Wrong label! "<<ID1<<"("<<PartiTags[ID1].first<<") "<<ID2<<"("<<PartiTags[ID2].first<<") "<<CHIndex[ID1][ID2].first<<"("<<dis<<") "<<CHIndex[ID1][ID2].second<<"("<<cnt<<")"<<endl;
//                    exit(1);
                }else if(CHIndex[ID1][ID2].second!=cnt ){
                    cout<<"Wrong cnt. "<<ID1<<"("<<PartiTags[ID1].first<<") "<<ID2<<"("<<PartiTags[ID2].first<<") "<<CHIndex[ID1][ID2].first<<"("<<dis<<") "<<CHIndex[ID1][ID2].second<<"("<<cnt<<")"<<endl;
                }
            }else{
                cout<<"Wrong. Not found. "<<ID1<<" "<<ID2<<" "<<dis<<endl; exit(1);
            }
        }
    }
}

//function of testing the throughput on random updates
void Graph::RandomUpdateThroughputTest(string updateFile, int batchNum, int batchSize, int batchInterval){
    bool ifDebug=false;
    ifDebug=true;
    int runtimes=10000;
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    cout<<"update file: "<<updateFile<<endl;
    string line;
    vector<string> vs;
    int ID1,ID2,oldW,newW,weight;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==1);
    int eNum=stoi(vs[0]);
    if(batchNum*batchSize>eNum){
        batchNum=eNum/batchSize;
        cout<<"Actual batch number: "<<batchNum<<endl;
    }
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNum);
    for(int i=0;i<batchNum;++i){
        for(int j=0;j<batchSize;++j){
            getline(IF,line);
            vs.clear();
            boost::split(vs,line,boost::is_any_of(" "));
            ID1=stoi(vs[0]), ID2=stoi(vs[1]), weight=stoi(vs[2]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    IF.close();
    cout<<"Update batch number: "<<batchNum<<" ; batch size: "<<batchSize<<" ; Batch interval: "<< batchInterval<<endl;

    string queryF = sourcePath+"/"+dataset + ".query";
    if(samePartiPortion!=-1){
        queryF=sourcePath+"/tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".sameParti_"+to_string(samePartiPortion)+".query";
    }
    ifstream IF2(queryF);
    if(!IF2){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF2>>num;
    for(int k=0;k<num;k++){
        IF2>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF2.close();

    //index maintenance
    Timer tt;
    Timer tRecord;
    double runT1=0, runT2 = 0;
    unsigned long long throughputNum=0;
    if(algoChoice==5){
        ProBeginVertexSetParti.assign(partiNum,vector<int>());
        vertexIDChLParti.assign(partiNum,set<int>());
        ProBeginVertexSetPartiExtend.assign(partiNum,vector<int>());
    }
    double updateTime=0;
    vector<double> stageUpdateT(5,0);
    vector<double> stageQueryT(5,0);
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;

    EffiStageCheck(ODpair, runtimes, stageQueryT);
//    EffiStageCheck(ODpair, 100, stageQueryT);
//    batchNum=2;
    overlayUpdateT=0;
    for(int i=0;i<batchNum;++i){
        stageDurations.assign(5,0);
        wBatchDec.clear(); wBatchInc.clear();
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[i].size();++j){
            ID1=batchUpdates[i][j].first.first, ID2=batchUpdates[i][j].first.second, weight=batchUpdates[i][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
//                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            if(j<batchSize/2){//decrease update
                for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                    if(it->first==ID2){
                        ifFind=true;
                        oldW=it->second;
//                        weight=(0.5+0.5*rand()/(RAND_MAX+1.0))*oldW;
//                        cout<<weight<<endl;
                        weight=0.5*oldW;
                        if(weight>0 && weight<oldW){
//                        cout<<"Dec "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                            wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                        }

                        break;
                    }
                }
            }
            else{//increase update
                for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                    if(it->first==ID2){
                        ifFind=true;
                        oldW=it->second;
//                        weight=(1+1*rand()/(RAND_MAX+1.0))*oldW;
                        weight=2*oldW;
                        if(weight>0 && weight>oldW) {
//                        cout<<"Inc "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                            wBatchInc.emplace_back(make_pair(ID1, ID2), make_pair(oldW, weight));
                        }
                        break;
                    }
                }
            }

            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),weight});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<weight<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<i<<" . Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;

        int id1=2, id2=3;
//        cout<<"Before update "<<i<<". "<<id1<<"("<<NodeOrder[id1]<<","<<PartiTags[id1].first<<") "<<id2<<"("<<NodeOrder[id2]<<","<<PartiTags[id2].first<<") "<<QueryPostMHL(id1,id2)<<"("<<Dijkstra(id1,id2,Neighbor)<<")"<<endl;
//        cout<<"Before update "<<i<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") "<<endl;
//        for(auto it=Tree[rank[id1]].vert.begin();it!=Tree[rank[id1]].vert.end();++it){
//            if(it->first==id2){
//                cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//                break;
//            }
//        }

        //Step 1: Decrease updates
        tt.start();
        if(!wBatchDec.empty()){
            cout<<"Decrease update. "<<wBatchDec.size()<<endl;
            tRecord.start();
//                    boost::thread_group thread;

            if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::DecBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT1)));
                DecBatchThroughputNP(wBatchDec, i, runT1);
            }else if(algoChoice==2){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
//                IndexMaintenance(wBatchDec, i, runT1);
            }

            else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                PMHLBatchUpdateDec(wBatchDec, i, runT1);
            }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                PostMHLBatchUpdateDec(wBatchDec, i, runT1);
            }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

//            if(ifDebug){
//                CorrectnessCheck(100);
//            }
        }

        //Step 2: Increase updates
        if(!wBatchInc.empty()){
            cout<<"Increase update. "<<wBatchInc.size()<<endl;
//            tRecord.start();
////            boost::thread_group thread;
            if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::IncBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT2)));
                IncBatchThroughputNP(wBatchInc, i, runT2);
            }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));
                PMHLBatchUpdateInc(wBatchInc, i, runT2);
            }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));//increase update is incorrect
                PostMHLBatchUpdateInc(wBatchInc, i, runT2);
            }

//            if(i==19){
//                ofstream OF(sourcePath+"/"+dataset+".CH19");
//                if(!OF.is_open()){
//                    cout<<"Fail to open file"<<endl; exit(1);
//                }
//                OF<<node_num<<" "<<edge_num<<endl;
//                for(int i=0;i<Neighbor.size();++i){
//                    for(auto it=Neighbor[i].begin();it!=Neighbor[i].end();++it){
//                        OF<<i<<" "<<it->first<<" "<<it->second<<endl;
//                    }
//                }
//                OF.close();
//                PostMHLIndexStoreCH(sourcePath+"/"+dataset+".CHIndex19");
//                cout<<"Write done."<<endl;
//            }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

//                    cout<<"After update: "<<endl;
//                    for(int i=0;i<Tree[rank[t]].disPost.size();++i){
//                        if(Tree[rank[t]].vAncestor[i]==hub){
//                            cout<<"Cnt: "<<Tree[rank[t]].disPost[i]<<" "<<Tree[rank[t]].cntPost[i]<<endl;
//                        }
//                    }
//                    for(auto it=Tree[rank[t]].vert.begin();it!=Tree[rank[t]].vert.end();++it){
//                        cout<<"sc: "<<t<<" "<<it->first<<" "<<it->second.first<<"("<<it->second.second<<")"<<endl;
//                    }
//                    cout<<"disB: "<<s<<" "<<bv<<" "<<Tree[rank[s]].disInf[BoundVertexMap[3][bv]]<<endl;


        }

        tt.stop();
        updateTime+=tt.GetRuntime();
        cout<<"Batch "<<i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;

        if(ifDebug){
            CorrectnessCheck(100);
        }

//        EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
        GetBatchThroughput(stageQueryT,batchInterval,throughputNum,stageUpdateT);
    }
//    for(int i=0;i<stageDurations.size();++i){
//        stageDurations[i]/=batchNum;
//    }
//    EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
    AverageStagePerformance(batchNum,stageUpdateT,stageQueryT);
    cout<<"\nPartiNum: "<<partiNum<<". Overall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average batch update Time: "<<updateTime/batchNum<<" s."<<endl;
}


//function of testing the throughput on random updates, for mix update, new
void Graph::RealUpdateThroughputTestQueueModel(string updateFile, int batchNum, double T_r, int workerNum, int regionNum) {//T_r is in seconds
    bool ifDebug=false;
//    ifDebug=true;
    int runtimes=10000;
//    runtimes=1000;
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    cout<<"Update file: "<<updateFile<<endl;
    string line;
    vector<string> vs;
    int ID1,ID2,oldW,newW;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==2);
    int batchNumMax=stoi(vs[0]);
    int batchInterval=stoi(vs[1]);
    int batchSize;
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNumMax);
    long long int aveBatchSize=0;
    int maxBatchSize=0, minBatchSize=INT32_MAX;
    for(int i=0;i<batchNumMax;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "));
        batchSize=stoi(vs[0]);
        assert(vs.size()==3*batchSize+1);
        aveBatchSize+=batchSize;
        if(maxBatchSize<batchSize) maxBatchSize=batchSize;
        if(minBatchSize>batchSize) minBatchSize=batchSize;
        for(int j=0;j<batchSize;++j){
            ID1=stoi(vs[3*j+1]), ID2=stoi(vs[3*j+2]), newW=stoi(vs[3*j+3]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),newW);
        }
    }
    IF.close();
    aveBatchSize/=batchNumMax;
    cout<<"Update batch number: "<<batchNumMax<<" ("<<batchNum<<") ; Average batch size: "<<aveBatchSize<<" ; Maximal batch size: "<<maxBatchSize<<" ; Minimal batch size: "<<minBatchSize<<" ; Batch interval: "<< batchInterval<<endl;
    vector<vector<pair<pair<int,int>,int>>> ODpairs;
    string fileQuery;
    fileQuery=sourcePath+"/"+dataset+".realQueries";
    ReadRealQueries(fileQuery, ODpairs);

    vector<pair<int,int>> ODpair;
    for(auto it=ODpairs[4].begin();it!=ODpairs[4].end();++it){
        ODpair.emplace_back(it->first.first,it->first.second);
    }
    cout<<"query number: "<<ODpair.size()<<endl;
//    vector<vector<int>> clusterP;
//    int clusterNum=4;
    if(algoChoice==6){
        PartitionRegionCluster(regionNum, clusterP, ODpair);
    }

    //index maintenance
    Timer tt;
    Timer tRecord;
    double runT1=0, runT2 = 0;
    unsigned long long throughputNum=0;
    if(algoChoice==5||algoChoice==6){
        ProBeginVertexSetParti.assign(partiNum,vector<int>());
        vertexIDChLParti.assign(partiNum,set<int>());
        ProBeginVertexSetPartiExtend.assign(partiNum,vector<int>());
        if(algoChoice==6){
            ProBeginVertexSetPartiInc.assign(partiNum,vector<int>());
            ProBeginVertexSetPartiExtendInc.assign(partiNum,vector<int>());
        }
    }
    double updateTime=0;
    vector<vector<double>> stageUpdateT(5,vector<double>());//(stage, update time)
    vector<vector<double>> stageQueryT(5,vector<double>());//(stage, query times)
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;
    vector<double> aveDuration(5,0.0);
    vector<double> aveQT(5,0);
//    EffiStageCheck(ODpair, runtimes, stageQueryT);//old
    EffiStageCheck(ODpair, ODpair.size(), stageQueryT, aveQT);
//    batchNum=2;
    overlayUpdateT=0;
    vector<unsigned long long> throughputNums;
    int query_i=0;
    int hourGap=4;
    int gap=3600*hourGap/batchInterval;
    lambdaPQueue.resize(1000001);
    lambdaCount.assign(1000001,0);
    int batchSizeSum=0;
    int count=0;
    double maxUTime=0;
    vector<int> sampleBatch;
    if(batchNum<=4 || batchNum==batchNumMax){
        for(int i=0;i<batchNum;++i){
            sampleBatch.push_back(i);
        }
    }else{
        sampleBatch=SampleBatches(sourcePath+"/"+dataset+"_20160105_"+to_string(batchInterval)+".batchUpdatesInfo",batchNum);
    }
    for(int i=0;i<sampleBatch.size();++i){
//    for(int batch_i=0;batch_i<batchNum;++batch_i){
//    for(int i=4;i<batchNum;++i){
//        if(batch_i%gap>1){
//        if(batch_i%gap>0){
//            continue;
//        }
        int batch_i=sampleBatch[i];
        stageDurations.assign(5,0);
        uStageDurations.assign(5,0);
        wBatchDec.clear(); wBatchInc.clear();
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[batch_i].size();++j){
            ID1=batchUpdates[batch_i][j].first.first, ID2=batchUpdates[batch_i][j].first.second, newW=batchUpdates[batch_i][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                if(it->first==ID2){
                    ifFind=true;
                    oldW=it->second;
//                        weight=max(oldW/2,1);
                    if(oldW>newW){
//                            if(i==4){
//                                cout<<"Dec "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<") "<<oldW<<" "<<weight<<endl;
//                            }
                        wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }else if(oldW<newW){
//                            if(i==4) {
//                                cout << "Inc " << ID1 << "(" << PartiTag[ID1].first << "," << PartiTag[ID1].second << ") " << ID2 << "(" << PartiTag[ID2].first << "," << PartiTag[ID2].second << ") " << oldW << " " << weight << endl;
//                            }
                        wBatchInc.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }
                    break;
                }
            }
            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),newW});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<newW<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<batch_i<<" . uNum: "<<wBatchDec.size()+wBatchInc.size()<<"("<<batchUpdates[batch_i].size()<<") ; Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;
        batchSizeSum=batchSizeSum+wBatchDec.size()+wBatchInc.size();
        count++;
        int id1=2, id2=3;


        updateBatchID=batch_i;
//            if(false){
        if(algoChoice==6){
            tt.start();
            cout<<"Mix update! Decrease size: "<<wBatchDec.size()<<" ; Increase size: "<<wBatchInc.size()<<endl;
//                VPLBatchUpdateMix(wBatchDec,wBatchInc,i,runT1);
            VPLBatchUpdateMixSchedule(wBatchDec,wBatchInc,batch_i,runT1);
            tt.stop();
        }
        else{
            //Step 1: Decrease updates
            tt.start();
            if(!wBatchDec.empty()){
                cout<<"Decrease update. "<<wBatchDec.size()<<endl;
                tRecord.start();
//                    boost::thread_group thread;

                if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::DecBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT1)));
                    DecBatchThroughputNP(wBatchDec, batch_i, runT1);
                }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                    PMHLBatchUpdateDec(wBatchDec, batch_i, runT1);
                }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                    PostMHLBatchUpdateDec(wBatchDec, batch_i, runT1);
                }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

//            if(ifDebug){
//                CorrectnessCheck(100);
//            }
            }

            //Step 2: Increase updates
            if(!wBatchInc.empty()){
                cout<<"Increase update. "<<wBatchInc.size()<<endl;
//            tRecord.start();
////            boost::thread_group thread;
                if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::IncBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT2)));
                    IncBatchThroughputNP(wBatchInc, batch_i, runT2);
                }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));
                    PMHLBatchUpdateInc(wBatchInc, batch_i, runT2);
                }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));
                    PostMHLBatchUpdateInc(wBatchInc, batch_i, runT2);
                }


            }

            tt.stop();
        }

        updateTime+=tt.GetRuntime();
        cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
        if(maxUTime<tt.GetRuntime()){
            maxUTime=tt.GetRuntime();
        }


        if(ifDebug){
            CorrectnessCheck(100);
        }
        for(int i=0;i<uStageDurations.size();++i){
            cout<<"Duration of U-Stage "<<i+1<<" : "<<uStageDurations[i]<<" s"<<endl;
            aveDuration[i]+=uStageDurations[i];
        }
//        for(int i=0;i<stageDurations.size()-1;++i){
//            cout<<"Duration of U-Stage "<<i+2<<" : "<<stageDurations[i]<<" s"<<endl;
//            aveDuration[i]+=stageDurations[i];
//        }


//        EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
        if(algoChoice==6){
            vector<double> duration = StageDurationCompute(batchInterval);
            double fastQT;
            for(int j=0;j<duration.size();++j){
//                stageQueryT[j].clear();//new
                stageUpdateT[j].push_back(duration[j]);
                if(duration[j]>0){
                    fastQT=aveQT[j];
                }
            }

            vector<pair<int,int>> ODpairThis;
            for(query_i=0;query_i<ODpairs[4].size();++query_i){
                if(ODpairs[4][query_i].second>=345600+batchInterval*batch_i){
                    if(ODpairs[4][query_i].second<345600+batchInterval*(batch_i+1)){
                        ODpairThis.emplace_back(ODpairs[4][query_i].first);
                    }else{
                        break;
                    }
                }
            }
            vector<vector<double>> stageQueryThis(5,vector<double>());//(stage, query times)
            EffiStageCheckReal(ODpairThis, ODpairThis.size(), stageQueryThis);
            for(int j=0;j<stageQueryThis.size();++j){
                for(auto it=stageQueryThis[j].begin();it!=stageQueryThis[j].end();++it){
                    stageQueryT[j].emplace_back(*it);
                }
            }

            vector<vector<double>> durations;
            StageDurationComputeSchedule(batchInterval,durations);
            throughputNum = ThroughputSimulateRealVPL(ODpairThis, stageQueryThis,durations,batchInterval,T_r,batchInterval/clusterP.size(),fastQT,workerNum);//2 batches update
            throughputNums.push_back(throughputNum);
        }
        else{//other algorithms
            vector<double> duration = StageDurationCompute(batchInterval);
            double fastQT;
            for(int j=0;j<duration.size();++j){
//                stageQueryT[j].clear();//new
                stageUpdateT[j].clear();
                stageUpdateT[j].push_back(duration[j]);
                if(duration[j]>0){
                    fastQT=aveQT[j];
                }
            }
            vector<pair<int,int>> ODpairThis;
            for(query_i=0;query_i<ODpairs[4].size();++query_i){
                if(ODpairs[4][query_i].second>=345600+batchInterval*batch_i){
                    if(ODpairs[4][query_i].second<345600+batchInterval*(batch_i+1)){
                        ODpairThis.emplace_back(ODpairs[4][query_i].first);
                    }else{
                        break;
                    }
                }
            }
            vector<vector<double>> stageQueryThis(5,vector<double>());//(stage, query times)
            EffiStageCheckReal(ODpairThis, ODpairThis.size(), stageQueryThis);
            for(int j=0;j<stageQueryThis.size();++j){
                for(auto it=stageQueryThis[j].begin();it!=stageQueryThis[j].end();++it){
                    stageQueryT[j].emplace_back(*it);
                }
            }

//            pair<double,double> resultE = ThroughputEstimate(stageQueryT, stageUpdateT, T_r, batchInterval);
//            cout<<"Estimated throughput: "<<resultE.first<<" ; fastest available query time: "<<resultE.second*1000<<" ("<<fastQT<<") ms"<<endl;
            throughputNum = ThroughputSimulateReal(stageQueryThis,stageUpdateT,batchInterval,T_r,batchInterval,fastQT,workerNum);//2 batches update
            throughputNums.push_back(throughputNum);
        }

    }

    cout<<"*************\nAverage batch update Time: "<<updateTime/batchNum<<" s."<<endl;
    for(int i=0;i<aveDuration.size();++i){
        cout<<"Average update time of U-Stage "<<i+1<<" : "<<aveDuration[i]/batchNum<<" s."<<endl;
    }
    for(int i=0;i<stageQueryT.size();++i){
        if(algoChoice==6 || algoChoice==5){//
            if(i<=PCH_No){
                cout<<"Average query time of Q-Stage "<<i+1<<" : "<<get_mean(stageQueryT[i])*1000<<" ms"<<endl;
            }else if(i>=PH2H_Post){
                cout<<"Average query time of Q-Stage "<<i<<" : "<<get_mean(stageQueryT[i])*1000<<" ms"<<endl;
            }
        }else{
            cout<<"Average query time of Q-Stage "<<i+1<<" : "<<get_mean(stageQueryT[i])*1000<<" ms"<<endl;
        }
    }

//    exit(0);
    if(algoChoice!=6){
        pair<double,double> resultE = ThroughputEstimate(stageQueryT, stageUpdateT, T_r, batchInterval);
        cout<<"Estimated throughput: "<<resultE.first<<" ; fastest available query time: "<<resultE.second*1000<<" ms"<<endl;
//        throughputNum = ThroughputSimulate(stageQueryT,stageUpdateT,batchInterval*2,T_r,batchInterval,resultE.second,workerNum);//2 batches update
//        cout<<"\nPartiNum: "<<partiNum<<". Throughput: "<<throughputNum<<" ; Average batch update Time: "<<updateTime/batchNum<<" s; average region update time: "<<updateTime/(batchNum*regionNum)<<endl;
        double aveThr=0;
        for(int i=0;i<throughputNums.size();++i){
            aveThr+=throughputNums[i];
        }
        cout<<"Batch number: "<<batchNum<<"("<<count<<"). Maximal update time: "<<maxUTime<<" ; average batch size: "<<batchSizeSum/count<<endl;
        cout<<"\nPartiNum: "<<partiNum<<". Average Throughput: "<<aveThr/throughputNums.size()<<" ; Average batch update Time: "<<updateTime/count<<" s"<<endl;
    }else{
        double aveThr=0;
        for(int i=0;i<throughputNums.size();++i){
            aveThr+=throughputNums[i];
        }
        cout<<"Average unaffected query number: "<<get_mean(unaffectedQNum)<<" ; total query number: "<<get_mean(totalQNum)<<endl;
        cout<<"Batch number: "<<batchNum<<"("<<count<<"). Maximal update time: "<<maxUTime<<" ; average batch size: "<<batchSizeSum/count<<endl;
        cout<<"\nPartiNum: "<<partiNum<<". Average Throughput: "<<aveThr/throughputNums.size()<<" ; Average batch update Time: "<<updateTime/count<<" s."<<endl;
    }


}

//function of testing the throughput on random updates
void Graph::RandomUpdateThroughputTestQueueModel(int batchNum, int batchSize, int batchInterval, double T_r, int workerNum) {//T_r is in seconds
    bool ifDebug=false;
//    ifDebug=true;
    int runtimes=10000;
//    runtimes=1000;
    string updateFile=sourcePath+"/"+dataset+".update";
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    cout<<"update file: "<<updateFile<<endl;
    string line;
    vector<string> vs;
    int ID1,ID2,oldW,newW,weight;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==1);
    int eNum=stoi(vs[0]);
    if(batchNum*batchSize>eNum){
        batchNum=eNum/batchSize;
        cout<<"Actual batch number: "<<batchNum<<endl;
    }
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNum);
    for(int i=0;i<batchNum;++i){
        for(int j=0;j<batchSize;++j){
            getline(IF,line);
            vs.clear();
            boost::split(vs,line,boost::is_any_of(" "));
            ID1=stoi(vs[0]), ID2=stoi(vs[1]), weight=stoi(vs[2]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    IF.close();
    cout<<"Update batch number: "<<batchNum<<" ; batch size: "<<batchSize<<" ; Batch interval: "<< batchInterval<<endl;

    string queryF = sourcePath+"/"+dataset + ".query";
    if(samePartiPortion!=-1){
        queryF=sourcePath+"/tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".sameParti_"+to_string(samePartiPortion)+".query";
    }
    ifstream IF2(queryF);
    if(!IF2){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF2>>num;
    for(int k=0;k<num;k++){
        IF2>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF2.close();

    //index maintenance
    Timer tt;
    Timer tRecord;
    double runT1=0, runT2 = 0;
    double throughputNum=0;
    if(algoChoice==5){
        ProBeginVertexSetParti.assign(partiNum,vector<int>());
        vertexIDChLParti.assign(partiNum,set<int>());
        ProBeginVertexSetPartiExtend.assign(partiNum,vector<int>());
    }
    double updateTime=0;
    vector<vector<double>> stageUpdateT(5,vector<double>());//(stage, update time)
    vector<vector<double>> stageQueryT(5,vector<double>());//(stage, query times)
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;
    vector<double> aveDuration(5,0.0);
    EffiStageCheck(ODpair, runtimes, stageQueryT);
//    EffiStageCheck(ODpair, runtimes, stageQueryT);
//    EffiStageCheck(ODpair, 100, stageQueryT);
//    batchNum=2;
    overlayUpdateT=0;
    for(int i=0;i<batchNum;++i){
        stageDurations.clear();
        stageDurations.assign(5,0);
        wBatchDec.clear(); wBatchInc.clear();
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[i].size();++j){
            ID1=batchUpdates[i][j].first.first, ID2=batchUpdates[i][j].first.second, weight=batchUpdates[i][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
//                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            if(j<batchSize/2){//decrease update
                for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                    if(it->first==ID2){
                        ifFind=true;
                        oldW=it->second;
//                        weight=(0.5+0.5*rand()/(RAND_MAX+1.0))*oldW;
//                        cout<<weight<<endl;
                        weight=oldW/2;
                        if(weight>0 && weight<oldW){
//                        cout<<"Dec "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                            wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                        }

                        break;
                    }
                }
            }
            else{//increase update
                for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                    if(it->first==ID2){
                        ifFind=true;
                        oldW=it->second;
//                        weight=(1+1*rand()/(RAND_MAX+1.0))*oldW;
                        weight=2*oldW;
                        if(weight>0 && weight>oldW) {
//                        cout<<"Inc "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                            wBatchInc.emplace_back(make_pair(ID1, ID2), make_pair(oldW, weight));
                        }
                        break;
                    }
                }
            }

            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),weight});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<weight<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<i<<" . Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;

        int id1=2, id2=3;
//        cout<<"Before update "<<i<<". "<<id1<<"("<<NodeOrder[id1]<<","<<PartiTags[id1].first<<") "<<id2<<"("<<NodeOrder[id2]<<","<<PartiTags[id2].first<<") "<<QueryPostMHL(id1,id2)<<"("<<Dijkstra(id1,id2,Neighbor)<<")"<<endl;
//        cout<<"Before update "<<i<<". "<<id1<<"("<<NodeOrder[id1]<<") "<<id2<<"("<<NodeOrder[id2]<<") "<<endl;
//        for(auto it=Tree[rank[id1]].vert.begin();it!=Tree[rank[id1]].vert.end();++it){
//            if(it->first==id2){
//                cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//                break;
//            }
//        }

        //Step 1: Decrease updates
        tt.start();
        if(!wBatchDec.empty()){
            cout<<"Decrease update. "<<wBatchDec.size()<<endl;
            tRecord.start();
//                    boost::thread_group thread;

            if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::DecBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT1)));
                DecBatchThroughputNP(wBatchDec, i, runT1);
            }else if(algoChoice==2){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
//                IndexMaintenance(wBatchDec, i, runT1);
            }

            else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                PMHLBatchUpdateDec(wBatchDec, i, runT1);
            }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                PostMHLBatchUpdateDec(wBatchDec, i, runT1);
            }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

//            if(ifDebug){
//                CorrectnessCheck(100);
//            }
        }

        //Step 2: Increase updates
        if(!wBatchInc.empty()){
            cout<<"Increase update. "<<wBatchInc.size()<<endl;
//            tRecord.start();
////            boost::thread_group thread;
            if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::IncBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT2)));
                IncBatchThroughputNP(wBatchInc, i, runT2);
            }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));
                PMHLBatchUpdateInc(wBatchInc, i, runT2);
            }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));//increase update is incorrect
                PostMHLBatchUpdateInc(wBatchInc, i, runT2);
            }

//            if(i==19){
//                ofstream OF(sourcePath+"/"+dataset+".CH19");
//                if(!OF.is_open()){
//                    cout<<"Fail to open file"<<endl; exit(1);
//                }
//                OF<<node_num<<" "<<edge_num<<endl;
//                for(int i=0;i<Neighbor.size();++i){
//                    for(auto it=Neighbor[i].begin();it!=Neighbor[i].end();++it){
//                        OF<<i<<" "<<it->first<<" "<<it->second<<endl;
//                    }
//                }
//                OF.close();
//                PostMHLIndexStoreCH(sourcePath+"/"+dataset+".CHIndex19");
//                cout<<"Write done."<<endl;
//            }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

//                    cout<<"After update: "<<endl;
//                    for(int i=0;i<Tree[rank[t]].disPost.size();++i){
//                        if(Tree[rank[t]].vAncestor[i]==hub){
//                            cout<<"Cnt: "<<Tree[rank[t]].disPost[i]<<" "<<Tree[rank[t]].cntPost[i]<<endl;
//                        }
//                    }
//                    for(auto it=Tree[rank[t]].vert.begin();it!=Tree[rank[t]].vert.end();++it){
//                        cout<<"sc: "<<t<<" "<<it->first<<" "<<it->second.first<<"("<<it->second.second<<")"<<endl;
//                    }
//                    cout<<"disB: "<<s<<" "<<bv<<" "<<Tree[rank[s]].disInf[BoundVertexMap[3][bv]]<<endl;


        }

        tt.stop();
        updateTime+=tt.GetRuntime();
        cout<<"Batch "<<i<<" . Update time: "<<tt.GetRuntime()<<" s."<<endl;

        if(ifDebug){
            CorrectnessCheck(100);
        }

        for(int i=0;i<stageDurations.size()-1;++i){
            cout<<"Duration of U-Stage "<<i+2<<" : "<<stageDurations[i]<<" s"<<endl;
            aveDuration[i]+=stageDurations[i];
        }


//        EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
        vector<double> duration = StageDurationCompute(batchInterval);
        for(int j=0;j<duration.size();++j){
            stageUpdateT[j].push_back(duration[j]);
        }

    }
//    for(int i=0;i<stageDurations.size();++i){
//        stageDurations[i]/=batchNum;
//    }
//    EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
//    AverageStagePerformance(batchNum,stageUpdateT,stageQueryT);
    cout<<"*************\nAverage batch update Time: "<<updateTime/batchNum<<" s."<<endl;
    for(int i=0;i<aveDuration.size()-1;++i){
        cout<<"Average update time of U-Stage "<<i+2<<" : "<<aveDuration[i]/batchNum<<" s ; query time: "<<get_mean(stageQueryT[i])*1000<<" ms"<<endl;
    }


    pair<double,double> resultE = ThroughputEstimate(stageQueryT, stageUpdateT, T_r, batchInterval);
    cout<<"Estimated throughput: "<<resultE.first<<" ; fastest available query time: "<<resultE.second*1000<<" ms"<<endl;
    throughputNum = ThroughputSimulate(stageQueryT,stageUpdateT,batchInterval*batchNum,T_r,batchInterval,resultE.second,workerNum);

    cout<<"\nPartiNum: "<<partiNum<<". Throughput: "<<throughputNum<<" ; Average batch update Time: "<<updateTime/batchNum<<" s."<<endl;
}

//function of testing the throughput of path-finding system, batchInterval is the time interval between two adjacent update batch (in seconds)
void Graph::SPThroughputTest(int updateType, bool ifBatch, int batchNum, int batchSize, int batchInterval, int runtimes) {
    cout<<"Shortest path query throughput test..."<<endl;
    // read updates
    string file = sourcePath+"/"+dataset + ".update";
    bool ifDebug=false;
    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand(0);
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);

    if(batchNum*batchSize>updateData.size()){
        batchNum=floor(updateData.size()/batchSize);
    }
    cout<<"Update batch: "<<batchNum<<" ; Batch size: "<<batchSize<<" ; Batch interval: "<< batchInterval<<endl;

    string queryF = sourcePath+"/"+dataset + ".query";
    if(samePartiPortion!=-1){
        queryF=sourcePath+"/tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".sameParti_"+to_string(samePartiPortion)+".query";
    }
    ifstream IF(queryF);
    if(!IF){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF.close();
    unsigned long long throughputNum=0;
    if(algoChoice==5){
        ProBeginVertexSetParti.assign(partiNum,vector<int>());
        vertexIDChLParti.assign(partiNum,set<int>());
        ProBeginVertexSetPartiExtend.assign(partiNum,vector<int>());
    }

    Timer tt;
    Timer tRecord;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 0:{
            break;
        }
        case 1:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;

            if(ifBatch){//for batch update
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
//                Graph g2=*this;
                vector<double> stageUpdateT(5,0);
                vector<double> stageQueryT(5,0);

                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    stageDurations.assign(5,0);
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[u*batchSize+i].first.first;
                        ID2 = updateData[u*batchSize+i].first.second;
                        oldW = updateData[u*batchSize+i].second;
                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
                        if(ifDebug){
                            //cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

//                    int s,t;
//                    s=112485, t=148515;//core-core
//                    s=41522, t=47616;//core-parti
//                    cout<<"Before update. "<<QueryPostMHLDebug(s,t)<<endl;

                    tRecord.start();
//                    boost::thread_group thread;

                    if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::DecBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT1)));
                        DecBatchThroughputNP(wBatch, u, runT1);
                    }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                        PMHLBatchUpdateDec(wBatch, u, runT1);
                    }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                        PostMHLBatchUpdateDec(wBatch, u, runT1);
                    }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

                    EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);

                    if(ifDebug){
//                        g2.CorrectnessCheckH2H(100);
                        CorrectnessCheck(100);
                    }
                }
                StagePerformanceShow(batchNum,stageUpdateT,stageQueryT);
                cout<<"\nPartiNum: "<<partiNum<<". Overall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average Decrease batch update Time: "<<runT1/batchNum<<" s."<<endl;
            }
            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            if(ifBatch){//for batch update
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                vector<double> stageUpdateT(5,0);
                vector<double> stageQueryT(5,0);

//                for(int u=5;u<6;++u){
//                for(int u=117;u<118;++u){
                for(int u=0;u<batchNum;++u){
                    wBatch.clear();
                    stageDurations.assign(5,0);
                    for(int i=0;i<batchSize;++i) {
                        ID1 = updateData[u * batchSize + i].first.first;
                        ID2 = updateData[u * batchSize + i].first.second;
                        oldW = updateData[u * batchSize + i].second;
                        newW = oldW * 1.5;
                        if (ifDebug) {
                            //cout << "Batch " << u << ": " << ID1 << " " << ID2 << " " << oldW << " " << newW << endl;
                        }
                        wBatch.emplace_back(make_pair(ID1, ID2), make_pair(oldW, newW));
                    }
//                    int s,t,hub;
//                    s=67432, t=67975;//same-parti, 13076
////                    s=115541, t=116033;//same-parti//115999
//                    cout<<"Before update. "<<QueryPostMHLDebug(s,t)<<endl;
//                    hub=s;
//                    for(int i=0;i<Tree[rank[t]].disPost.size();++i){
//                        if(Tree[rank[t]].vAncestor[i]==hub){
//                            cout<<"Cnt: "<<Tree[rank[t]].disPost[i]<<" "<<Tree[rank[t]].cntPost[i]<<endl;
//                        }
//                    }
//                    for(auto it=Tree[rank[t]].vert.begin();it!=Tree[rank[t]].vert.end();++it){
//                        cout<<"sc: "<<t<<" "<<it->first<<" "<<it->second.first<<"("<<it->second.second<<")"<<endl;
//                    }
//                    int bv=11922;
//                    cout<<"disB: "<<s<<" "<<bv<<" "<<Tree[rank[s]].disInf[BoundVertexMap[3][bv]]<<endl;

                    tRecord.start();
                    boost::thread_group thread;
                    if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::IncBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT2)));
                        IncBatchThroughputNP(wBatch, u, runT2);
                    }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));
                        PMHLBatchUpdateInc(wBatch, u, runT2);
                    }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));//increase update is incorrect
                        PostMHLBatchUpdateInc(wBatch, u, runT2);
                    }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

                    EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);

//                    cout<<"After update: "<<endl;
//                    for(int i=0;i<Tree[rank[t]].disPost.size();++i){
//                        if(Tree[rank[t]].vAncestor[i]==hub){
//                            cout<<"Cnt: "<<Tree[rank[t]].disPost[i]<<" "<<Tree[rank[t]].cntPost[i]<<endl;
//                        }
//                    }
//                    for(auto it=Tree[rank[t]].vert.begin();it!=Tree[rank[t]].vert.end();++it){
//                        cout<<"sc: "<<t<<" "<<it->first<<" "<<it->second.first<<"("<<it->second.second<<")"<<endl;
//                    }
//                    cout<<"disB: "<<s<<" "<<bv<<" "<<Tree[rank[s]].disInf[BoundVertexMap[3][bv]]<<endl;

                    if(ifDebug){
                        CorrectnessCheck(100);
                    }
//                    exit(0);
                }
                StagePerformanceShow(batchNum,stageUpdateT,stageQueryT);
                cout<<"\nPartiNum: "<<partiNum<<". Overall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average Increase batch update Time: "<<runT2/batchNum<<" s."<<endl;
            }


            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

void Graph::StagePerformanceShow(int batchNum, vector<double>& stageUpdateT, vector<double>& stageQueryT){
    if(algoChoice==1){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]/batchNum<<" ms."<<endl;
        cout<<"Stage 2 (CH). Average update time: "<<stageUpdateT[CH]/batchNum<<" s; average query time: "<<1000*stageQueryT[CH]/batchNum<<" ms."<<endl;
        cout<<"Stage 3 (H2H). Average duration: "<<stageUpdateT[H2H]/batchNum<<" s; average query time: "<<1000*stageQueryT[H2H]/batchNum<<" ms."<<endl;
    }else if(algoChoice==3){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]/batchNum<<" ms."<<endl;
        cout<<"Stage 2 (PCH). Average update time: "<<stageUpdateT[PCH_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PCH_No]/batchNum<<" ms."<<endl;
        cout<<"Stage 3 (No-boundary PMHL). Average update time: "<<stageUpdateT[PH2H_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_No]/batchNum<<" ms."<<endl;
        cout<<"Stage 4 (Post-boundary PMHL). Average update time: "<<stageUpdateT[PH2H_Post]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Post]/batchNum<<" ms."<<endl;
        cout<<"Stage 5 (Cross-boundary PMHL). Average duration: "<<stageUpdateT[PH2H_Cross]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Cross]/batchNum<<" ms."<<endl;
    }else if(algoChoice==5){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]/batchNum<<" ms."<<endl;
        cout<<"Stage 2 (PCH). Average update time: "<<stageUpdateT[PCH_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PCH_No]/batchNum<<" ms."<<endl;
        cout<<"Stage 3 (Post-boundary PostMHL). Average update time: "<<stageUpdateT[PH2H_Post]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Post]/batchNum<<" ms."<<endl;
        cout<<"Stage 4 (Cross-boundary PostMHL). Average duration: "<<stageUpdateT[PH2H_Cross]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Cross]/batchNum<<" ms."<<endl;
    }
}

void Graph::AverageStagePerformance(int batchNum, vector<double> &stageUpdateT, vector<double> &stageQueryT) {
    if(algoChoice==1){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]<<" ms."<<endl;
        cout<<"Stage 2 (CH). Average update time: "<<stageUpdateT[CH]/batchNum<<" s; average query time: "<<1000*stageQueryT[CH]<<" ms."<<endl;
        cout<<"Stage 3 (H2H). Average duration: "<<stageUpdateT[H2H]/batchNum<<" s; average query time: "<<1000*stageQueryT[H2H]<<" ms."<<endl;
    }else if(algoChoice==3){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]<<" ms."<<endl;
        cout<<"Stage 2 (PCH). Average update time: "<<stageUpdateT[PCH_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PCH_No]<<" ms."<<endl;
        cout<<"Stage 3 (No-boundary PMHL). Average update time: "<<stageUpdateT[PH2H_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_No]<<" ms."<<endl;
        cout<<"Stage 4 (Post-boundary PMHL). Average update time: "<<stageUpdateT[PH2H_Post]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Post]<<" ms."<<endl;
        cout<<"Stage 5 (Cross-boundary PMHL). Average duration: "<<stageUpdateT[PH2H_Cross]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Cross]<<" ms."<<endl;
    }else if(algoChoice==5){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]<<" ms."<<endl;
        cout<<"Stage 2 (PCH). Average update time: "<<stageUpdateT[PCH_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PCH_No]<<" ms."<<endl;
        cout<<"Stage 3 (Overlay label update). Average duration: "<<overlayUpdateT/batchNum<<" s."<<endl;
        cout<<"Stage 4 (Post-boundary PostMHL). Average update time: "<<stageUpdateT[PH2H_Post]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Post]<<" ms."<<endl;
        cout<<"Stage 5 (Cross-boundary PostMHL). Average duration: "<<stageUpdateT[PH2H_Cross]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Cross]<<" ms."<<endl;
    }
}

//function for throughput test of decrease update
void Graph::PMHLBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
    if(!partiBatch.empty()){
        if(threadnum==1){
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
                DecreasePartiBatchUpdateCheckCH(pid, it->second, overlayBatch, false, updatedSCs[pid]);//old
            }
        }
        else{
            if(partiBatch.size()>threadnum){
                cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
                vector<vector<int>> processID;
                processID.assign(threadnum, vector<int>());
                vector<int> vertices;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    vertices.emplace_back(pid);
                }
                ThreadDistribute(vertices, processID);
                boost::thread_group thread;
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheckCHV, this, processID[j], boost::ref(partiBatch), boost::ref(overlayBatch), false, boost::ref(updatedSCs) ));
                }
                thread.join_all();
            }else{
                boost::thread_group thread;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false, boost::ref(updatedSCs[pid]) ));//old
//            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), true, boost::ref(updatedSCs[pid]) ));
                }
                thread.join_all();
            }

        }

    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }

            if(newdis<olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                sm->notify();
            }
        }
    }
//    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<<updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);//weightOverlay collect the changed edges on overlay graph
    }

    cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,false);
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]+=tt2.GetRuntime();
    tt2.start();
//    cout<<"algoQuery: PCH-No"<<endl;

    if(threadnum==1){
        DecreaseOverlayBatchLabel(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL
        if(!partiBatch.empty()){
//        if(partiBatch.size()>threadnum){
//            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
//        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
//            cout<<"partition "<<pid<<" "<<heightMaxs[pid]<< endl;
                DecreasePartiBatchLabel(Trees[pid],ranks[pid], heightMaxs[pid], ProBeginVertexSetParti[pid], vertexIDChLParti[pid]);
            }
        }

    }else{

        if(partiBatch.size()+1>threadnum){
            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::DecreaseOverlayBatchLabel, this, boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(ProBeginVertexSetOverlay), boost::ref(vertexIDChLOverlay) ));

            vector<vector<int>> processID;
            processID.assign(threadnum-1, vector<int>());
            vector<int> vertices;
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);

            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchLabelV, this, processID[j], boost::ref(Trees), boost::ref(ranks), boost::ref(heightMaxs), boost::ref(ProBeginVertexSetParti), boost::ref(vertexIDChLParti) ));
            }
            thread.join_all();

        }else{
            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::DecreaseOverlayBatchLabel, this, boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(ProBeginVertexSetOverlay), boost::ref(vertexIDChLOverlay) ));

            if(!partiBatch.empty()){
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
//            cout<<"partition "<<pid<<" "<<heightMaxs[pid]<< endl;
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchLabel, this, boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], boost::ref(ProBeginVertexSetParti[pid]), boost::ref(vertexIDChLParti[pid]) ));
                }
            }
            thread.join_all();
        }

    }


    algoQuery=PH2H_No;
    tt2.stop();
    stageDurations[PCH_No]+=tt2.GetRuntime();
    // repair the partition index
    if(algoUpdate>=PH2H_Post){
//        cout<<"repair post-boundary index"<<endl;
        tt2.start();
        if(threadnum==1){
            Repair_PartiIndex(false, false, partiBatch);//post
        }else{
            Repair_PartiIndex(true, false, partiBatch);//post
        }
        algoQuery=PH2H_Post;
        tt2.stop();
        stageDurations[PH2H_No]+=tt2.GetRuntime();
        if(algoUpdate==PH2H_Cross){
            tt2.start();
//            RefreshExtensionLabels(partiBatch);//extend
            if(threadnum==1){
                RefreshExtensionLabelsNoAllPair(partiBatch, false);//extend
            }else{
                RefreshExtensionLabelsNoAllPair(partiBatch,true);//extend
            }

            algoQuery=PH2H_Cross;
            tt2.stop();
            stageDurations[PH2H_Post]+=tt2.GetRuntime();
            tt2.start();
        }
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

//function for throughput test of decrease update
void Graph::PMHLBatchUpdateDecOpt(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
//            cout<<"Parti update: "<<pid1<<". "<<a<<" "<<b<<" "<<newW<<endl;
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
    if(!partiBatch.empty()){
//        if(partiBatch.size()>threadnum){
//            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
//        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), true , boost::ref(updatedSCs[pid])));
        }
        thread.join_all();
    }

    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }

            if(newdis<olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                sm->notify();
            }
        }
    }
    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<<updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);//weightOverlay collect the changed edges on overlay graph
    }
    cout<<"overlayBatch size: "<<overlayBatch.size()<<endl;
    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,false);
    algoQuery=PCH_No;
//    cout<<"algoQuery: PCH-No"<<endl;
    if(algoUpdate>PCH_No){
        DecreaseOverlayBatchLabel(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL
    }
    // repair the partition index
    if(algoUpdate>=PH2H_Post){

        Repair_PartiIndex(true, false, partiBatch);//post
//        Repair_PartiIndexForOpt(false, false, partiBatch);//post


        algoQuery=PH2H_Post;
        if(algoUpdate==PH2H_Cross){
//            RefreshExtensionLabels(partiBatch);//extend
            RefreshExtensionLabelsNoAllPair(partiBatch,true);//extend
            algoQuery=PH2H_Cross;
        }
    }

    tt.stop();
    runT1+=tt.GetRuntime();
    cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
}

//function for VPL's throughput test of mix update
/*void Graph::VPLBatchUpdateMixSchedule(vector<pair<pair<int,int>,pair<int,int>>>& wBatchDec, vector<pair<pair<int,int>,pair<int,int>>>& wBatchInc, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    Timer tt3;
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchDec; partiBatchDec.clear();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchInc; partiBatchInc.clear();
    map<int,vector<pair<pair<int, int>, pair<int, int>>>> overlayBatchDecM;
    map<int,vector<pair<pair<int, int>, pair<int, int>>>> overlayBatchIncM;
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchDec;
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchInc;
//    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchDecPure;
//    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchIncPure;

    unordered_set<int> scheduleP;
    for(int k=0;k<wBatchDec.size();k++) {
        int a = wBatchDec[k].first.first;
        int b = wBatchDec[k].first.second;
        int oldW = wBatchDec[k].second.first;
        int newW = wBatchDec[k].second.second;

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(PartiTag[a].second && PartiTag[b].second){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatchDec.emplace_back(wBatchDec[k]);
//            overlayBatchDec.insert({wBatchDec[k].first,wBatchDec[k].second});
//            scheduleP.insert(partiNum);
        }else{
            if(!PartiTag[a].second && PartiTag[b].second){//b is boundary
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchDec[pid1].emplace_back(wBatchDec[k]);
                scheduleP.insert(pid1);
            }else if(PartiTag[a].second && !PartiTag[b].second){//a is boundary
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchDec[pid2].emplace_back(wBatchDec[k]);
                scheduleP.insert(pid2);
            }
            else{
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    partiBatchDec[pid1].emplace_back(wBatchDec[k]);
                    scheduleP.insert(pid1);
                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }
            }
        }
    }
    for(int k=0;k<wBatchInc.size();k++) {
        int a = wBatchInc[k].first.first;
        int b = wBatchInc[k].first.second;
        int oldW = wBatchInc[k].second.first;
        int newW = wBatchInc[k].second.second;

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(PartiTag[a].second && PartiTag[b].second){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatchInc.emplace_back(wBatchInc[k]);
//            overlayBatchInc.insert({wBatchInc[k].first,wBatchInc[k].second});
//            scheduleP.insert(partiNum);
        }else{
            if(!PartiTag[a].second && PartiTag[b].second){//b is boundary
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchInc[pid1].emplace_back(wBatchInc[k]);
                scheduleP.insert(pid1);
            }else if(PartiTag[a].second && !PartiTag[b].second){//a is boundary
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchInc[pid2].emplace_back(wBatchInc[k]);
                scheduleP.insert(pid2);
            }
            else{
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    partiBatchInc[pid1].emplace_back(wBatchInc[k]);
                    scheduleP.insert(pid1);
                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }

            }
        }
    }
    cout<<"Partition number for update: "<<scheduleP.size()<<endl;
    bool ifExt=false;
    for(auto it=overlayBatchDec.begin();it!=overlayBatchDec.end();++it){
        int pa = PartiTag[it->first.first].first;
        int pb = PartiTag[it->first.second].first;
        if(scheduleP.find(pa)!=scheduleP.end()){//if found
            overlayBatchDecM[pa].emplace_back(it->first,it->second);
        }else if(scheduleP.find(pb)!=scheduleP.end()) {//if found
            overlayBatchDecM[pb].emplace_back(it->first, it->second);
        }else{
            overlayBatchDecM[pa].emplace_back(it->first, it->second);
            scheduleP.insert(pa);
            ifExt=true;
        }
    }
    for(auto it=overlayBatchInc.begin();it!=overlayBatchInc.end();++it){
        int pa = PartiTag[it->first.first].first;
        int pb = PartiTag[it->first.second].first;
        if(scheduleP.find(pa)!=scheduleP.end()){//if found
            overlayBatchIncM[pa].emplace_back(it->first,it->second);
        }else if(scheduleP.find(pb)!=scheduleP.end()) {//if found
            overlayBatchIncM[pb].emplace_back(it->first, it->second);
        }else{
            overlayBatchIncM[pa].emplace_back(it->first, it->second);
            scheduleP.insert(pa);
            ifExt=true;
        }
    }
    if(ifExt){
        cout<<"Extended Partition number for update: "<<scheduleP.size()<<endl;
    }
//    vector<vector<int>> clusterP;
//    PartitionUpdateClustering(scheduleP, clusterP);
    bool ifFirst=true;
    int scheduleNum=0;
    for(auto it=scheduleP.begin();it!=scheduleP.end();++it){
        cout<<"Schedule "<<scheduleNum<<endl;
        scheduleNum++;
        tt3.start();
        algoQuery=Dijk;
        tt2.start();
        map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchDecSchedule;
        map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchIncSchedule;
        vector<pair<pair<int, int>, pair<int, int>>> overlayBatchDecSchedule;
        vector<pair<pair<int, int>, pair<int, int>>> overlayBatchIncSchedule;
        if(overlayBatchDecM.find(*it)!=overlayBatchDecM.end()){
            overlayBatchDecSchedule=overlayBatchDecM[*it];
            for(int k=0;k<overlayBatchDecSchedule.size();++k){
                int a = overlayBatchDecSchedule[k].first.first;
                int b = overlayBatchDecSchedule[k].first.second;
                int oldW = overlayBatchDecSchedule[k].second.first;
                int newW = overlayBatchDecSchedule[k].second.second;

                for(int i=0;i<Neighbor[a].size();i++){
                    if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                        assert(oldW==Neighbor[a][i].second);
                        if(newW<Neighbor[a][i].second){
                            Neighbor[a][i].second=newW;
                        }else{
                            cout<<"Wrong update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                        }
                        break;
                    }
                }
                for(int i=0;i<Neighbor[b].size();i++){
                    if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                        assert(oldW==Neighbor[b][i].second);
                        if(newW<Neighbor[b][i].second){
                            Neighbor[b][i].second=newW;
                        }else{
                            cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                        }

                        break;
                    }
                }
            }
        }
        if(overlayBatchIncM.find(*it)!=overlayBatchIncM.end()){
            overlayBatchIncSchedule=overlayBatchIncM[*it];
            for(int k=0;k<overlayBatchIncSchedule.size();++k){
                int a = overlayBatchIncSchedule[k].first.first;
                int b = overlayBatchIncSchedule[k].first.second;
                int oldW = overlayBatchIncSchedule[k].second.first;
                int newW = overlayBatchIncSchedule[k].second.second;
                for(int i=0;i<Neighbor[a].size();i++){
                    if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                        assert(oldW==Neighbor[a][i].second);
                        if(newW>Neighbor[a][i].second){
                            Neighbor[a][i].second=newW;
                        }else{
                            cout<<"Wrong update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                        }
                        break;
                    }
                }
                for(int i=0;i<Neighbor[b].size();i++){
                    if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                        assert(oldW==Neighbor[b][i].second);
                        if(newW>Neighbor[b][i].second){
                            Neighbor[b][i].second=newW;
                        }else{
                            cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                        }
                        break;
                    }
                }
            }
        }
        if(partiBatchDec.find(*it)!=partiBatchDec.end()){//if found
            partiBatchDecSchedule.insert({*it,partiBatchDec[*it]});
            for(int k=0;k<partiBatchDec[*it].size();++k){
                int a = partiBatchDec[*it][k].first.first;
                int b = partiBatchDec[*it][k].first.second;
                int oldW = partiBatchDec[*it][k].second.first;
                int newW = partiBatchDec[*it][k].second.second;
                for(int i=0;i<Neighbor[a].size();i++){
                    if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                        assert(oldW==Neighbor[a][i].second);
                        if(newW<Neighbor[a][i].second){
                            Neighbor[a][i].second=newW;
                        }else{
                            cout<<"Wrong update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                        }
                        break;
                    }
                }
                for(int i=0;i<Neighbor[b].size();i++){
                    if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                        assert(oldW==Neighbor[b][i].second);
                        if(newW<Neighbor[b][i].second){
                            Neighbor[b][i].second=newW;
                        }else{
                            cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                        }

                        break;
                    }
                }
            }
        }
        if(partiBatchInc.find(*it)!=partiBatchInc.end()){//if found
            partiBatchIncSchedule.insert({*it,partiBatchInc[*it]});
            for(int k=0;k<partiBatchInc[*it].size();++k){
                int a = partiBatchInc[*it][k].first.first;
                int b = partiBatchInc[*it][k].first.second;
                int oldW = partiBatchInc[*it][k].second.first;
                int newW = partiBatchInc[*it][k].second.second;
                for(int i=0;i<Neighbor[a].size();i++){
                    if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                        assert(oldW==Neighbor[a][i].second);
                        if(newW>Neighbor[a][i].second){
                            Neighbor[a][i].second=newW;
                        }else{
                            cout<<"Wrong update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                        }
                        break;
                    }
                }
                for(int i=0;i<Neighbor[b].size();i++){
                    if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                        assert(oldW==Neighbor[b][i].second);
                        if(newW>Neighbor[b][i].second){
                            Neighbor[b][i].second=newW;
                        }else{
                            cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                        }
                        break;
                    }
                }
            }
        }

        innerAffectedParti.assign(partiNum,true);
        outerAffectedParti.assign(partiNum,false);
        AffectedParti.assign(partiNum, true);
        affectedParti.clear(); affectedPartiInc.clear();
        overlayShortcutDec.clear();
        vUpdated.assign(node_num, false);
        ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();
        ProBeginVertexSetOverlayInc.clear();
        for(int i=0;i<partiNum;++i){
            ProBeginVertexSetParti[i].clear(); vertexIDChLParti[i].clear();
            ProBeginVertexSetPartiInc[i].clear();
            ProBeginVertexSetPartiExtend[i].clear();
        }
        for(int i=0;i<Tree.size();++i){
            Tree[i].DisRe.clear();
            Tree[i].DisRePost.clear();
        }

        double tDetect=0;
        double tSCU=0;
        boost::thread_group thread1;
        thread1.add_thread(new boost::thread(&Graph::VPLInnerUnaffectedPartiDetectByOracle, this, boost::ref(partiBatchDecSchedule), boost::ref(partiBatchIncSchedule),boost::ref(overlayBatchDecSchedule),boost::ref(overlayBatchIncSchedule),boost::ref(tDetect)));
        thread1.add_thread(new boost::thread(&Graph::VPLShortcutUpdate, this, boost::ref(partiBatchDecSchedule), boost::ref(partiBatchIncSchedule),boost::ref(overlayBatchDecSchedule),boost::ref(overlayBatchIncSchedule),boost::ref(tSCU)));
        thread1.join_all();
        if(tDetect>tSCU){
            cout<<"Shortcut update is faster than oracle detection!!! "<<tDetect<<" "<<tSCU<<endl;
        }

//        VPLInnerUnaffectedPartiDetectByOracle(partiBatchDecSchedule, partiBatchIncSchedule, overlayBatchDecSchedule, overlayBatchIncSchedule);
//        // Step 1: Update shortcuts
//        // Approach 1: partitioned version
//        tt.start();
//        /// Decrease update the shortcuts of non-boundary vertices
//        DecreasePartiBatchUpdateCheckVPL(partiBatchDecSchedule, overlayBatchDecSchedule, true);//multi-thread
////    DecreasePartiBatchUpdateCheckVPL(partiBatch, overlayBatch, false);//single-thread
////        cout<<"Decrease overlayBatch size: "<<overlayBatchDecSchedule.size()<<endl;
//        /// Inrease update the shortcuts of non-boundary vertices
//        IncreasePartiBatchUpdateCheckVPL(partiBatchIncSchedule, overlayBatchIncSchedule, true);//multi-thread, bottom-up shortcut update
////    IncreasePartiBatchUpdateCheckVPL(partiBatchInc, overlayBatchInc, false);//multi-thread, bottom-up shortcut update
////        cout<<"Increase overlayBatch size: "<<overlayBatchIncSchedule.size()<<endl;
//        tt.stop();
//        cout<<"Partition shortcut update time: "<<tt.GetRuntime()<<" s."<<endl;
//
//        tt.start();
////    VPLOverlayGraphUpdate(overlayBatchDec,overlayBatchInc,false);
//        VPLOverlayGraphUpdate(overlayBatchDecSchedule,overlayBatchIncSchedule,true);
//        tt.stop();
//        cout<<"Overlay graph update time: "<<tt.GetRuntime()<<" s."<<endl;
////    CorrectnessCheckCore(100);
//
//        tt.start();
//        /// Overlay shortcut update
//        DecreaseOverlayBatchVPL(overlayBatchDecSchedule,Tree,rank,heightMax,false);
//        IncreaseOverlayBatchVPL(overlayBatchIncSchedule,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
//        tt.stop();
//        cout<<"Overlay shortcut update time: "<<tt.GetRuntime()<<" s."<<endl;


//        CorrectnessCheck(100);

        // Approach 2: non-partitioned version
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,true);//direct bottom-up, with label construction
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,false);//direct bottom-up, without label construction
        algoQuery=PCH_No;
        tt2.stop();
        stageDurations[Dijk]+=tt2.GetRuntime();
//    cout<<"algoQuery: PCH-No"<<endl;

        //parallel version
        double tOverlay=0,tPost=0,tCross=0;
        if(algoUpdate==PH2H_Post) {
            tt2.start();
            DecreaseOverlayBatchLabelVPL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlayInc, VidtoTNid);
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            overlayUpdateT+=tOverlay;
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            tt2.start();
            // Step 3: Update partition index
            if(!partiBatchDecSchedule.empty() ||!partiBatchIncSchedule.empty()) {
                tt.start();
                BoundShortcutsCheck(true, 3);// check which partition is affected and update the boundary shortcuts
                tt.stop();
                cout << "All-pair boundary check time: " << tt.GetRuntime() << " s. Decrease affectedParti size: "
                     << affectedParti.size() <<" ; Increase affectedParti size: "<<affectedPartiInc.size()<< endl;
                ifRepaired.assign(partiNum, false);
            }


            /// Identify the unaffected partitions
            double tDetect=0;
            VPLUnaffectedPartiDetectTruth(overlayBatchDecSchedule,overlayBatchIncSchedule,tDetect);

            MixRepair_PartiIndexVPLPost(true, partiBatchDecSchedule, partiBatchIncSchedule, tPost);//decrease
            tt2.stop();
            stageDurations[PCH_No] = stageDurations[PCH_No]+ tt2.GetRuntime() + tOverlay;
        }
        else if(algoUpdate==PH2H_Cross){
            tt2.start();
            DecreaseOverlayBatchLabelVPL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlayInc, VidtoTNid);
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            overlayUpdateT+=tOverlay;

            // Step 3: Update partition index
            if(!partiBatchDecSchedule.empty() ||!partiBatchIncSchedule.empty()) {
                tt.start();
                BoundShortcutsCheck(true, 3);// check which partition is affected and update the boundary shortcuts
                tt.stop();
                cout << "All-pair boundary check time: " << tt.GetRuntime() << " s. Decrease affectedParti size: "
                     << affectedParti.size() <<" ; Increase affectedParti size: "<<affectedPartiInc.size()<< endl;
                ifRepaired.assign(partiNum, false);
            }
            /// Identify the unaffected partitions
            double tDetect=0;
            VPLUnaffectedPartiDetectTruth(overlayBatchDecSchedule,overlayBatchIncSchedule,tDetect);
//            if(true){
            if(threadnum<=affectedParti.size() && threadnum<=affectedPartiInc.size()) {
//                cout<<"No enough thread for parallelization between post-boundary and cross-boundary."<<endl;
                MixRepair_PartiIndexVPLPost(true, partiBatchDecSchedule, partiBatchIncSchedule, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, false, partiBatch);//post
                tt2.stop();
                stageDurations[PCH_No] = tt2.GetRuntime();
                cout<<"Post-boundary update time: "<<tPost<<" s"<<endl;
                if (algoUpdate == PH2H_Cross) {
                    tt2.start();
                    MixRefreshCrossLabelsVPL(true, tCross, threadnum);
//                RefreshCrossLabelsVPL(true, false, tCross, threadnum);//extend, parallel
//            RefreshExtensionLabelsPostMHL(false, false);//extend, single-thread
                    tt2.stop();
                    cout<<"Cross-boundary update time: "<<tCross<<" s"<<endl;
                    stageDurations[PH2H_Post] = tt2.GetRuntime();
                }
            }else{
//                cout<<"Post-boundary and cross-boundary could be parallelized!!!"<<endl;
                boost::thread_group thread;
                thread.add_thread(new boost::thread(&Graph::MixRepair_PartiIndexVPLPost, this, true, boost::ref(partiBatchDecSchedule), boost::ref(partiBatchIncSchedule), boost::ref(tPost) ));
//                thread.add_thread(new boost::thread(&Graph::MixRefreshCrossLabelsVPL, this, true, boost::ref(tCross), threadnum-max(affectedParti.size(),affectedPartiInc.size())));
                thread.add_thread(new boost::thread(&Graph::MixRefreshCrossLabelsVPL, this, true, boost::ref(tCross), threadnum-1));
                thread.join_all();
                if(tPost<tCross){
                    stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tPost;
                    stageDurations[PH2H_Post]=stageDurations[PH2H_Post]+tCross-tPost;
                    cout<<"Post and Cross-boundary update time: "<<tPost<<" "<<tCross<<endl;
                }else{
                    cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
                    stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tCross;
//            stageDurations[PH2H_Post]=0;
                }

            }

        }
        tt3.stop();
        CorrectnessCheck(10);
//        cout<<"Update time for partition "<<*it<<" : "<<tt3.GetRuntime()<<" s."<<endl;

    }

//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}*/


//function for VPL's throughput test of mix update, with cluster
void Graph::VPLBatchUpdateMixSchedule(vector<pair<pair<int,int>,pair<int,int>>>& wBatchDec, vector<pair<pair<int,int>,pair<int,int>>>& wBatchInc, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    Timer tt3;
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchDec; partiBatchDec.clear();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchInc; partiBatchInc.clear();
    map<int,vector<pair<pair<int, int>, pair<int, int>>>> overlayBatchDecM;
    map<int,vector<pair<pair<int, int>, pair<int, int>>>> overlayBatchIncM;
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchDec;
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchInc;
//    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchDecPure;
//    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchIncPure;

    unordered_set<int> scheduleP;
    for(int k=0;k<wBatchDec.size();k++) {
        int a = wBatchDec[k].first.first;
        int b = wBatchDec[k].first.second;
        int oldW = wBatchDec[k].second.first;
        int newW = wBatchDec[k].second.second;

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(PartiTag[a].second && PartiTag[b].second){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatchDec.emplace_back(wBatchDec[k]);
//            overlayBatchDec.insert({wBatchDec[k].first,wBatchDec[k].second});
//            scheduleP.insert(partiNum);
        }else{
            if(!PartiTag[a].second && PartiTag[b].second){//b is boundary
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchDec[pid1].emplace_back(wBatchDec[k]);
                scheduleP.insert(pid1);
            }else if(PartiTag[a].second && !PartiTag[b].second){//a is boundary
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchDec[pid2].emplace_back(wBatchDec[k]);
                scheduleP.insert(pid2);
            }
            else{
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    partiBatchDec[pid1].emplace_back(wBatchDec[k]);
                    scheduleP.insert(pid1);
                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }
            }
        }
    }
    for(int k=0;k<wBatchInc.size();k++) {
        int a = wBatchInc[k].first.first;
        int b = wBatchInc[k].first.second;
        int oldW = wBatchInc[k].second.first;
        int newW = wBatchInc[k].second.second;

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(PartiTag[a].second && PartiTag[b].second){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatchInc.emplace_back(wBatchInc[k]);
//            overlayBatchInc.insert({wBatchInc[k].first,wBatchInc[k].second});
//            scheduleP.insert(partiNum);
        }else{
            if(!PartiTag[a].second && PartiTag[b].second){//b is boundary
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchInc[pid1].emplace_back(wBatchInc[k]);
                scheduleP.insert(pid1);
            }else if(PartiTag[a].second && !PartiTag[b].second){//a is boundary
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchInc[pid2].emplace_back(wBatchInc[k]);
                scheduleP.insert(pid2);
            }
            else{
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    partiBatchInc[pid1].emplace_back(wBatchInc[k]);
                    scheduleP.insert(pid1);
                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }

            }
        }
    }
    cout<<"Partition number for update: "<<scheduleP.size()<<" ("<< partiNum<<")"<< endl;
    bool ifExt=false;
    for(auto it=overlayBatchDec.begin();it!=overlayBatchDec.end();++it){
        int pa = PartiTag[it->first.first].first;
        int pb = PartiTag[it->first.second].first;
        if(scheduleP.find(pa)!=scheduleP.end()){//if found
            overlayBatchDecM[pa].emplace_back(it->first,it->second);
        }else if(scheduleP.find(pb)!=scheduleP.end()) {//if found
            overlayBatchDecM[pb].emplace_back(it->first, it->second);
        }else{
            overlayBatchDecM[pa].emplace_back(it->first, it->second);
            scheduleP.insert(pa);
            ifExt=true;
        }
    }
    for(auto it=overlayBatchInc.begin();it!=overlayBatchInc.end();++it){
        int pa = PartiTag[it->first.first].first;
        int pb = PartiTag[it->first.second].first;
//        if(it->first.first==131800 && it->first.second==132139){
//            cout<<"Find here. "<<it->first.first<<" "<<it->first.second<<" "<<it->second.first<<" "<<it->second.second<<endl;
//        }
        if(scheduleP.find(pa)!=scheduleP.end()){//if found
            overlayBatchIncM[pa].emplace_back(it->first,it->second);
        }else if(scheduleP.find(pb)!=scheduleP.end()) {//if found
            overlayBatchIncM[pb].emplace_back(it->first, it->second);
        }else{
            overlayBatchIncM[pa].emplace_back(it->first, it->second);
            scheduleP.insert(pa);
            ifExt=true;
        }
    }
    if(ifExt){
        cout<<"Extended Partition number for update: "<<scheduleP.size()<<endl;
    }
//    vector<vector<int>> clusterP;
//    PartitionUpdateClustering(scheduleP, clusterP);
//    vector<int> temp;
//    for(auto it=scheduleP.begin();it!=scheduleP.end();++it){
//        temp.push_back(*it);
//    }
//    clusterP.push_back(temp);
    bool ifFirst=true;
    int scheduleNum=0;
    VPLUpdateTimes.assign(clusterP.size(),vector<double>(5,0));
    partiAffectInfo.assign(clusterP.size(),vector<pair<bool,int>>(partiNum,make_pair(true,0)));
    for(int ci=0;ci<clusterP.size();++ci){
//    for(auto it=scheduleP.begin();it!=scheduleP.end();++it){
        cout<<"Schedule "<<scheduleNum<<endl;
        scheduleNum++;
        tt3.start();

        tt2.start();
        map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchDecSchedule;
        map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchIncSchedule;
        vector<pair<pair<int, int>, pair<int, int>>> overlayBatchDecSchedule;
        vector<pair<pair<int, int>, pair<int, int>>> overlayBatchIncSchedule;
        for(auto it=clusterP[ci].begin();it!=clusterP[ci].end();++it) {
            if (overlayBatchDecM.find(*it) != overlayBatchDecM.end()) {
                for (int k = 0; k < overlayBatchDecM[*it].size(); ++k) {
                    int a = overlayBatchDecM[*it][k].first.first;
                    int b = overlayBatchDecM[*it][k].first.second;
                    int oldW = overlayBatchDecM[*it][k].second.first;
                    int newW = overlayBatchDecM[*it][k].second.second;
                    overlayBatchDecSchedule.emplace_back(overlayBatchDecM[*it][k].first,overlayBatchDecM[*it][k].second);
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                            assert(oldW == Neighbor[a][i].second);
                            if (newW < Neighbor[a][i].second) {
                                Neighbor[a][i].second = newW;
                            } else {
                                cout << "Wrong update. " << a << " " << b << " " << Neighbor[a][i].second << " " << newW
                                     << endl;
                                exit(1);
                            }
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                            assert(oldW == Neighbor[b][i].second);
                            if (newW < Neighbor[b][i].second) {
                                Neighbor[b][i].second = newW;
                            } else {
                                cout << "Wrong decrease update. " << a << " " << b << " " << Neighbor[a][i].second
                                     << " " << newW << endl;
                                exit(1);
                            }

                            break;
                        }
                    }
                }
            }
            if (overlayBatchIncM.find(*it) != overlayBatchIncM.end()) {
                for (int k = 0; k < overlayBatchIncM[*it].size(); ++k) {
                    int a = overlayBatchIncM[*it][k].first.first;
                    int b = overlayBatchIncM[*it][k].first.second;
                    int oldW = overlayBatchIncM[*it][k].second.first;
                    int newW = overlayBatchIncM[*it][k].second.second;
                    overlayBatchIncSchedule.emplace_back(overlayBatchIncM[*it][k].first, overlayBatchIncM[*it][k].second);
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                            assert(oldW == Neighbor[a][i].second);
                            if (newW > Neighbor[a][i].second) {
                                Neighbor[a][i].second = newW;
                            } else {
                                cout << "Wrong update. " << a << " " << b << " " << Neighbor[a][i].second << " " << newW
                                     << endl;
                                exit(1);
                            }
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                            assert(oldW == Neighbor[b][i].second);
                            if (newW > Neighbor[b][i].second) {
                                Neighbor[b][i].second = newW;
                            } else {
                                cout << "Wrong decrease update. " << a << " " << b << " " << Neighbor[a][i].second
                                     << " " << newW << endl;
                                exit(1);
                            }
                            break;
                        }
                    }
                }
            }
            if (partiBatchDec.find(*it) != partiBatchDec.end()) {//if found
                partiBatchDecSchedule.insert({*it, partiBatchDec[*it]});
                for (int k = 0; k < partiBatchDec[*it].size(); ++k) {
                    int a = partiBatchDec[*it][k].first.first;
                    int b = partiBatchDec[*it][k].first.second;
                    int oldW = partiBatchDec[*it][k].second.first;
                    int newW = partiBatchDec[*it][k].second.second;
//                    if((a==99549 && b==81831)||(a==81831 && b==99549)){
//                        cout<<"Find. "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl;
//                    }
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                            assert(oldW == Neighbor[a][i].second);
                            if (newW < Neighbor[a][i].second) {
                                Neighbor[a][i].second = newW;
                            } else {
                                cout << "Wrong update. " << a << " " << b << " " << Neighbor[a][i].second << " " << newW
                                     << endl;
                                exit(1);
                            }
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                            assert(oldW == Neighbor[b][i].second);
                            if (newW < Neighbor[b][i].second) {
                                Neighbor[b][i].second = newW;
                            } else {
                                cout << "Wrong decrease update. " << a << " " << b << " " << Neighbor[a][i].second
                                     << " " << newW << endl;
                                exit(1);
                            }

                            break;
                        }
                    }
                }
            }
            if (partiBatchInc.find(*it) != partiBatchInc.end()) {//if found
                partiBatchIncSchedule.insert({*it, partiBatchInc[*it]});
                for (int k = 0; k < partiBatchInc[*it].size(); ++k) {
                    int a = partiBatchInc[*it][k].first.first;
                    int b = partiBatchInc[*it][k].first.second;
                    int oldW = partiBatchInc[*it][k].second.first;
                    int newW = partiBatchInc[*it][k].second.second;
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                            assert(oldW == Neighbor[a][i].second);
                            if (newW > Neighbor[a][i].second) {
                                Neighbor[a][i].second = newW;
                            } else {
                                cout << "Wrong update. " << a << " " << b << " " << Neighbor[a][i].second << " " << newW
                                     << endl;
                                exit(1);
                            }
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                            assert(oldW == Neighbor[b][i].second);
                            if (newW > Neighbor[b][i].second) {
                                Neighbor[b][i].second = newW;
                            } else {
                                cout << "Wrong decrease update. " << a << " " << b << " " << Neighbor[a][i].second
                                     << " " << newW << endl;
                                exit(1);
                            }
                            break;
                        }
                    }
                }
            }
        }
        tt2.stop();
        algoQuery=Dijk;
        uStageDurations[Dijk]+=tt2.GetRuntime();
        cout<<"partiDec size: "<<partiBatchDecSchedule.size()<<" ; partiInc size: "<<partiBatchIncSchedule.size()<<" ; overlayDec size: "<<overlayBatchDecSchedule.size()<<" ; overlayInc size: "<<overlayBatchIncSchedule.size()<<endl;
        tt2.start();
        innerAffectedParti.assign(partiNum,true);
        outerAffectedParti.assign(partiNum,false);
        AffectedParti.assign(partiNum, true);
        affectedParti.clear(); affectedPartiInc.clear();
        overlayShortcutDec.clear();
        vUpdated.assign(node_num, false);
        ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();
        ProBeginVertexSetOverlayInc.clear();
        overlayLDec.clear(), overlayLInc.clear();
        partiLDec.assign(partiNum,map<pair<int,int>,pair<int,int>>());
        partiLInc.assign(partiNum,map<pair<int,int>,int>());
        checkNumParti.assign(partiNum,0);
        for(int i=0;i<partiNum;++i){
            ProBeginVertexSetParti[i].clear(); vertexIDChLParti[i].clear();
            ProBeginVertexSetPartiInc[i].clear();
            ProBeginVertexSetPartiExtend[i].clear();
        }
        for(int i=0;i<Tree.size();++i){
            Tree[i].DisRe.clear();
            Tree[i].DisRePost.clear();
            Tree[i].disChange.assign(Tree[i].dis.size(),0);
            Tree[i].ifU=false;
//            Tree[i].disOldM.clear();
        }

        double tDetect=0;
        double tSCU=0;
        int uNum=0;
//        VPLInnerUnaffectedPartiDetectByOracle(partiBatchDecSchedule, partiBatchIncSchedule,overlayBatchDecSchedule,overlayBatchIncSchedule,uNum,tDetect,partiAffectInfo[ci]);
//        VPLShortcutUpdate(partiBatchDecSchedule, partiBatchIncSchedule,overlayBatchDecSchedule,overlayBatchIncSchedule,tSCU);
        boost::thread_group thread1;
        thread1.add_thread(new boost::thread(&Graph::VPLInnerUnaffectedPartiDetectByOracle, this, boost::ref(partiBatchDecSchedule), boost::ref(partiBatchIncSchedule),boost::ref(overlayBatchDecSchedule),boost::ref(overlayBatchIncSchedule), boost::ref(uNum), boost::ref(tDetect), boost::ref(partiAffectInfo[ci])));
        thread1.add_thread(new boost::thread(&Graph::VPLShortcutUpdate, this, boost::ref(partiBatchDecSchedule), boost::ref(partiBatchIncSchedule),boost::ref(overlayBatchDecSchedule),boost::ref(overlayBatchIncSchedule),boost::ref(tSCU)));
        thread1.join_all();
        cout<<"Shortcut update time: "<<tSCU<<" s ; Oracle-based detect time: "<<tDetect<<" s ; Inner-unaffected partition number by oracle: "<<uNum<<endl;
        tt2.stop();
        algoQuery=PCH_No;
        uStageDurations[PCH_No]+=tt2.GetRuntime();
        stageDurations[Dijk]+=tt2.GetRuntime();
        VPLUpdateTimes[ci][Dijk]=tt2.GetRuntime();
//    cout<<"algoQuery: PCH-No"<<endl;

        //parallel version
        double tOverlay=0,tPost=0,tCross=0;
        if(algoUpdate>PCH_No){
            tt2.start();
//            DecreaseOverlayBatchLabelVPL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
//            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlayInc, VidtoTNid);
            MixVPLOverlayLabelUpdateTopDown(overlayLDec,overlayLInc,Tree,rank,VidtoTNid);//correct version by top-down mechanism
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            overlayUpdateT+=tOverlay;
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            uStageDurations[PH2H_No] += tOverlay;

            if(algoUpdate==PH2H_Post) {
                tt2.start();
                // Step 3: Update partition index
                if(!partiBatchDecSchedule.empty() ||!partiBatchIncSchedule.empty()) {
                    tt.start();
                    BoundShortcutsCheck(true, MixUpdate);// check which partition is affected and update the boundary shortcuts
                    tt.stop();
                    cout << "All-pair boundary check time: " << tt.GetRuntime() << " s. Decrease affectedParti size: "
                         << affectedParti.size() <<" ; Increase affectedParti size: "<<affectedPartiInc.size()<< endl;
                    ifRepaired.assign(partiNum, false);
                }


                /// Identify the unaffected partitions
                double tDetect=0;
                VPLUnaffectedPartiDetectTruth(overlayBatchDecSchedule,overlayBatchIncSchedule,tDetect, partiAffectInfo[ci]);
                MixRepair_PartiIndexVPLPost(true, partiBatchDecSchedule, partiBatchIncSchedule, tPost);//decrease
                tt2.stop();
                uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tt2.GetRuntime();
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tt2.GetRuntime();
                VPLUpdateTimes[ci][PCH_No]=tOverlay+tt2.GetRuntime();
            }
            else if(algoUpdate==PH2H_Cross){
                // Step 3: Update partition index
                double tBound=0;
                if(!partiBatchDecSchedule.empty() ||!partiBatchIncSchedule.empty()) {
                    tt.start();
                    BoundShortcutsCheck(true, MixUpdate);// check which partition is affected and update the boundary shortcuts
                    tt.stop();
                    tBound=tt.GetRuntime();
                    cout << "All-pair boundary check time: " << tt.GetRuntime() << " s. Decrease affectedParti size: "
                         << affectedParti.size() <<" ; Increase affectedParti size: "<<affectedPartiInc.size()<< endl;
                    ifRepaired.assign(partiNum, false);
                }
                /// Identify the unaffected partitions
                double tDetect=0;
                VPLUnaffectedPartiDetectTruth(overlayBatchDecSchedule,overlayBatchIncSchedule,tDetect,partiAffectInfo[ci]);
//            if(true){
                if(threadnum<=affectedParti.size() && threadnum<=affectedPartiInc.size()) {
//                cout<<"No enough thread for parallelization between post-boundary and cross-boundary."<<endl;
                    tt2.start();
                    MixRepair_PartiIndexVPLPost(true, partiBatchDecSchedule, partiBatchIncSchedule, tPost);//post
                    tt2.stop();
                    uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tt2.GetRuntime()+tBound;
                    stageDurations[PCH_No] = stageDurations[PCH_No]+tt2.GetRuntime()+tBound+tOverlay;
                    VPLUpdateTimes[ci][PCH_No]=tt2.GetRuntime()+tBound+tOverlay;
                    cout<<"Post-boundary update time: "<<tPost<<" s"<<endl;
                    if (algoUpdate == PH2H_Cross) {
                        tt2.start();
                        MixRefreshCrossLabelsVPL(true, tCross, threadnum);
                        tt2.stop();
                        cout<<"Cross-boundary update time: "<<tCross<<" s"<<endl;
                        uStageDurations[PH2H_Cross] += tt2.GetRuntime();
                        stageDurations[PH2H_Post] += tt2.GetRuntime();
                    }
                }
                else{
//                cout<<"Post-boundary and cross-boundary could be parallelized!!!"<<endl;
                    boost::thread_group thread;
                    thread.add_thread(new boost::thread(&Graph::MixRepair_PartiIndexVPLPost, this, true, boost::ref(partiBatchDecSchedule), boost::ref(partiBatchIncSchedule), boost::ref(tPost) ));
//                thread.add_thread(new boost::thread(&Graph::MixRefreshCrossLabelsVPL, this, true, boost::ref(tCross), threadnum-max(affectedParti.size(),affectedPartiInc.size())));
                    thread.add_thread(new boost::thread(&Graph::MixRefreshCrossLabelsVPL, this, true, boost::ref(tCross), threadnum-2));
                    thread.join_all();
                    if(tPost<tCross){
                        uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tPost+tBound;
                        uStageDurations[PH2H_Cross]=uStageDurations[PH2H_Cross]+tCross-tPost;
                        stageDurations[PCH_No] = stageDurations[PCH_No]+tPost+tBound+tOverlay;
                        VPLUpdateTimes[ci][PCH_No]=tPost+tBound+tOverlay;
                        stageDurations[PH2H_Post]=stageDurations[PH2H_Post]+tCross-tPost;
                        VPLUpdateTimes[ci][PH2H_Post]=tCross-tPost;
                        cout<<"Post and Cross-boundary update time: "<<tPost+tBound<<" "<<tCross<<endl;
                    }else{
                        cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
                        uStageDurations[PH2H_Cross] = uStageDurations[PH2H_Cross]+tCross+tBound;
                        stageDurations[PCH_No] = stageDurations[PCH_No]+tCross+tBound+ tOverlay;
                        VPLUpdateTimes[ci][PCH_No]=tCross+tBound+ tOverlay;
//            stageDurations[PH2H_Post]=0;
                    }

                }

            }
        }

        tt3.stop();
//        CorrectnessCheck(10);
//        cout<<"Update time for partition "<<*it<<" : "<<tt3.GetRuntime()<<" s."<<endl;

    }

//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

void Graph::VPLShortcutUpdate(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatchDecSchedule, map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatchIncSchedule, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchDecSchedule, vector<pair<pair<int, int>, pair<int, int>>>& overlayBatchIncSchedule, double& tSc){
    Timer tt;
    Timer tt2;
    tt2.start();
    // Step 1: Update shortcuts
    // Approach 1: partitioned version
//    tt.start();
//    DecreasePartiBatchUpdateCheckVPL(partiBatchDecSchedule, overlayBatchDecSchedule, true);/// Decrease update the shortcuts of non-boundary vertices
//    IncreasePartiBatchUpdateCheckVPL(partiBatchIncSchedule, overlayBatchIncSchedule, true);/// Increase update the shortcuts of non-boundary vertices
//    tt.stop();
//    cout<<"Partition shortcut update time: "<<tt.GetRuntime()<<" s."<<endl;
    MixVPLBatchShortcutUpdateParti(partiBatchDecSchedule,partiBatchIncSchedule,overlayBatchDecSchedule,overlayBatchIncSchedule,true);///Mix update
//    cout<<"overlayBatchDec size: "<<overlayBatchDecSchedule.size()<<" ; overlayBatchInc size: "<<overlayBatchIncSchedule.size()<<endl;
//    tt.start();
//    VPLOverlayGraphUpdate(overlayBatchDecSchedule,overlayBatchIncSchedule,true);///new remove
//    tt.stop();
//    if(tt.GetRuntime()>0.001){
//        cout<<"Overlay graph update time: "<<tt.GetRuntime()<<" s."<<endl;
//    }
//    CorrectnessCheckCore(100);
    /// Overlay shortcut update
//    tt.start();
//    DecreaseOverlayBatchVPL(overlayBatchDecSchedule,Tree,rank,heightMax,false);
//    IncreaseOverlayBatchVPL(overlayBatchIncSchedule,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    MixVPLBatchShortcutUpdateOverlay(overlayBatchDecSchedule,overlayBatchIncSchedule,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
//    tt.stop();
//    cout<<"Overlay shortcut update time: "<<tt.GetRuntime()<<" s."<<endl;
    tt2.stop();
    tSc=tt2.GetRuntime();
//    cout<<"Shortcut update time: "<<tSc<<" s."<<endl;
}

//function for VPL's throughput test of mix update
void Graph::VPLBatchUpdateMix(vector<pair<pair<int,int>,pair<int,int>>>& wBatchDec, vector<pair<pair<int,int>,pair<int,int>>>& wBatchInc, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchDec; partiBatchDec.clear();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatchInc; partiBatchInc.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchDec;
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatchInc;
    innerAffectedParti.assign(partiNum,true);
    outerAffectedParti.assign(partiNum,false);
    AffectedParti.assign(partiNum, true);
    affectedParti.clear(); affectedPartiInc.clear();
    overlayShortcutDec.clear();
    for(int k=0;k<wBatchDec.size();k++) {
        int a = wBatchDec[k].first.first;
        int b = wBatchDec[k].first.second;
        int oldW = wBatchDec[k].second.first;
        int newW = wBatchDec[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                assert(oldW==Neighbor[a][i].second);
                if(newW<Neighbor[a][i].second){
                    Neighbor[a][i].second=newW;
                }else{
                    cout<<"Wrong update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                assert(oldW==Neighbor[b][i].second);
                if(newW<Neighbor[b][i].second){
                    Neighbor[b][i].second=newW;
                }else{
                    cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }

                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(PartiTag[a].second && PartiTag[b].second){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatchDec.emplace_back(wBatchDec[k]);
        }else{
            if(!PartiTag[a].second && PartiTag[b].second){//b is boundary
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchDec[pid1].emplace_back(wBatchDec[k]);

            }else if(PartiTag[a].second && !PartiTag[b].second){//a is boundary
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchDec[pid2].emplace_back(wBatchDec[k]);
            }
            else{
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    partiBatchDec[pid1].emplace_back(wBatchDec[k]);
                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }
            }
        }
    }
    for(int k=0;k<wBatchInc.size();k++) {
        int a = wBatchInc[k].first.first;
        int b = wBatchInc[k].first.second;
        int oldW = wBatchInc[k].second.first;
        int newW = wBatchInc[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                assert(oldW==Neighbor[a][i].second);
                if(newW>Neighbor[a][i].second){
                    Neighbor[a][i].second=newW;
                }else{
                    cout<<"Wrong update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                assert(oldW==Neighbor[b][i].second);
                if(newW>Neighbor[b][i].second){
                    Neighbor[b][i].second=newW;
                }else{
                    cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(PartiTag[a].second && PartiTag[b].second){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatchInc.emplace_back(wBatchInc[k]);
        }else{
            if(!PartiTag[a].second && PartiTag[b].second){//b is boundary
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchInc[pid1].emplace_back(wBatchInc[k]);
            }else if(PartiTag[a].second && !PartiTag[b].second){//a is boundary
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                partiBatchInc[pid2].emplace_back(wBatchInc[k]);
            }
            else{
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    partiBatchInc[pid1].emplace_back(wBatchInc[k]);

                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }

            }
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();
    ProBeginVertexSetOverlayInc.clear();
    for(int i=0;i<partiNum;++i){
        ProBeginVertexSetParti[i].clear(); vertexIDChLParti[i].clear();
        ProBeginVertexSetPartiInc[i].clear();
        ProBeginVertexSetPartiExtend[i].clear();
    }
    for(int i=0;i<Tree.size();++i){
        Tree[i].DisRe.clear();
        Tree[i].DisRePost.clear();
    }
    // Step 1: Update shortcuts
    // Approach 1: partitioned version
    tt.start();
    /// Decrease update the shortcuts of non-boundary vertices
    DecreasePartiBatchUpdateCheckVPL(partiBatchDec, overlayBatchDec, true);//multi-thread
//    DecreasePartiBatchUpdateCheckVPL(partiBatch, overlayBatch, false);//single-thread
    cout<<"Decrease overlayBatch size: "<<overlayBatchDec.size()<<endl;
    /// Inrease update the shortcuts of non-boundary vertices
    IncreasePartiBatchUpdateCheckVPL(partiBatchInc, overlayBatchInc, true);//multi-thread, bottom-up shortcut update
//    IncreasePartiBatchUpdateCheckVPL(partiBatchInc, overlayBatchInc, false);//multi-thread, bottom-up shortcut update
    cout<<"Increase overlayBatch size: "<<overlayBatchInc.size()<<endl;
    tt.stop();
    cout<<"Partition shortcut update time: "<<tt.GetRuntime()<<" s."<<endl;


    tt.start();
//    VPLOverlayGraphUpdate(overlayBatchDec,overlayBatchInc,false);
    VPLOverlayGraphUpdate(overlayBatchDec,overlayBatchInc,true);
    tt.stop();
    cout<<"Overlay graph update time: "<<tt.GetRuntime()<<" s."<<endl;
//    CorrectnessCheckCore(100);


    tt.start();
    /// Overlay shortcut update
    DecreaseOverlayBatchVPL(overlayBatchDec,Tree,rank,heightMax,false);
    IncreaseOverlayBatchVPL(overlayBatchInc,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    tt.stop();
    cout<<"Overlay shortcut update time: "<<tt.GetRuntime()<<" s."<<endl;


//    CorrectnessCheck(100);

    // Approach 2: non-partitioned version
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,true);//direct bottom-up, with label construction
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,false);//direct bottom-up, without label construction
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]+=tt2.GetRuntime();
//    cout<<"algoQuery: PCH-No"<<endl;
    if(true){
//    if(threadnum<=affectedParti.size() && threadnum<=affectedPartiInc.size()){
        cout<<"No enough thread for parallelization between post-boundary and cross-boundary."<<endl;
        //unparalleled version
        if(algoUpdate>=PH2H_Post){
            double tOverlay=0,tPost=0,tCross=0;
            tt2.start();
            // Step 2: Update overlay index
            DecreaseOverlayBatchLabelVPL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlayInc, VidtoTNid);
            tt2.stop();
            tOverlay=tt2.GetRuntime();
            cout<<"Overlay index update: "<<tOverlay<<endl;

            tt2.start();
            // Step 3: Update partition index
            if(!partiBatchDec.empty() ||!partiBatchInc.empty()) {
                tt.start();
                BoundShortcutsCheck(true, MixUpdate);// check which partition is affected and update the boundary shortcuts
                tt.stop();
                cout << "All-pair boundary check time: " << tt.GetRuntime() << " s. Decrease affectedParti size: "
                     << affectedParti.size() <<" ; Increase affectedParti size: "<<affectedPartiInc.size()<< endl;
                ifRepaired.assign(partiNum, false);
            }
            /// Identify the inner-unaffected partitions
            tt.start();
            VPLInnerUnaffectedPartiDetect(overlayBatchDec, overlayBatchInc, false);
            tt.stop();
            double tDetect=tt.GetRuntime();
            cout<<"Inner-unaffected partition detect time: "<<tDetect<<" s."<<endl;
            /// Identify the outer-unaffected partitions
            tt2.start();
            VPLOuterUnaffectedPartiDetect();
            tt2.stop();
            cout<<"Outer-unaffected partition identify time: "<<tt2.GetRuntime()<<" s."<<endl;
//            CorrectnessCheck(100);
//            tt.start();
////    BPFlagUpdateIncrease(overlayBatchInc, true);
//////    BPFlagUpdateIncrease(overlayBatchInc, false);
////    tt.stop();
////    cout<<"Increase BP-Flag update time: "<<tt.GetRuntime()<<" s."<<endl;
////    tt.start();
////    BPFlagUpdateDecrease(overlayBatchDec, true);
//////    BPFlagUpdateDecrease(overlayBatchDec, false);
//            VPLBPFlagConstruct(true);
//            tt.stop();
//            cout<<"Decrease BP-Flag update time: "<<tt.GetRuntime()<<" s."<<endl;

            MixRepair_PartiIndexVPLPost(true, partiBatchDec, partiBatchInc, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, false, partiBatch);//post
            tt2.stop();
            stageDurations[PCH_No]=tt2.GetRuntime();
            if(algoUpdate==PH2H_Cross){
                tt2.start();
                MixRefreshCrossLabelsVPL(true,tCross,threadnum);
//                RefreshCrossLabelsVPL(true, false, tCross, threadnum);//extend, parallel
//            RefreshExtensionLabelsPostMHL(false, false);//extend, single-thread
                tt2.stop();
                stageDurations[PH2H_Post]=tt2.GetRuntime();
            }
        }
    }
    else{//if thread number is sufficient
        cout<<"Post-boundary and cross-boundary could be parallelized!!!"<<endl;
        //parallel version
        double tOverlay=0,tPost=0,tCross=0;
        if(algoUpdate==PH2H_Post) {
            tt2.start();
            DecreaseOverlayBatchLabelVPL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlayInc, VidtoTNid);
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            overlayUpdateT+=tOverlay;
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            tt2.start();
            // Step 3: Update partition index
            if(!partiBatchDec.empty() ||!partiBatchInc.empty()) {
                tt.start();
                BoundShortcutsCheck(true, MixUpdate);// check which partition is affected and update the boundary shortcuts
                tt.stop();
                cout << "All-pair boundary check time: " << tt.GetRuntime() << " s. Decrease affectedParti size: "
                     << affectedParti.size() <<" ; Increase affectedParti size: "<<affectedPartiInc.size()<< endl;
                ifRepaired.assign(partiNum, false);
            }
            /// Identify the inner-unaffected partitions
            tt.start();
            VPLInnerUnaffectedPartiDetect(overlayBatchDec, overlayBatchInc, false);
            tt.stop();
            double tDetect=tt.GetRuntime();
            cout<<"Inner-unaffected partition detect time: "<<tDetect<<" s."<<endl;
            /// Identify the outer-unaffected partitions
            tt2.start();
            VPLOuterUnaffectedPartiDetect();
            tt2.stop();
            cout<<"Outer-unaffected partition identify time: "<<tt2.GetRuntime()<<" s."<<endl;
//            CorrectnessCheck(100);
//            tt.start();
////    BPFlagUpdateIncrease(overlayBatchInc, true);
//////    BPFlagUpdateIncrease(overlayBatchInc, false);
////    tt.stop();
////    cout<<"Increase BP-Flag update time: "<<tt.GetRuntime()<<" s."<<endl;
////    tt.start();
////    BPFlagUpdateDecrease(overlayBatchDec, true);
//////    BPFlagUpdateDecrease(overlayBatchDec, false);
//            VPLBPFlagConstruct(true);
//            tt.stop();
//            cout<<"Decrease BP-Flag update time: "<<tt.GetRuntime()<<" s."<<endl;

            MixRepair_PartiIndexVPLPost(true, partiBatchDec, partiBatchInc, tPost);//decrease
            tt2.stop();
            stageDurations[PCH_No] = stageDurations[PCH_No]+ tt2.GetRuntime() + tOverlay;
        }else if(algoUpdate==PH2H_Cross){
            tt2.start();
            DecreaseOverlayBatchLabelVPL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlayInc, VidtoTNid);
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            overlayUpdateT+=tOverlay;

            // Step 3: Update partition index
            if(!partiBatchDec.empty() ||!partiBatchInc.empty()) {
                tt.start();
                BoundShortcutsCheck(true, MixUpdate);// check which partition is affected and update the boundary shortcuts
                tt.stop();
                cout << "All-pair boundary check time: " << tt.GetRuntime() << " s. Decrease affectedParti size: "
                     << affectedParti.size() <<" ; Increase affectedParti size: "<<affectedPartiInc.size()<< endl;
                ifRepaired.assign(partiNum, false);
            }
            /// Identify the inner-unaffected partitions
            tt.start();
            VPLInnerUnaffectedPartiDetect(overlayBatchDec, overlayBatchInc, false);
            tt.stop();
            double tDetect=tt.GetRuntime();
            cout<<"Inner-unaffected partition detect time: "<<tDetect<<" s."<<endl;
            /// Identify the outer-unaffected partitions
            tt2.start();
            VPLOuterUnaffectedPartiDetect();
            tt2.stop();
            cout<<"Outer-unaffected partition identify time: "<<tt2.GetRuntime()<<" s."<<endl;
//            CorrectnessCheck(100);
//            tt.start();
////    BPFlagUpdateIncrease(overlayBatchInc, true);
//////    BPFlagUpdateIncrease(overlayBatchInc, false);
////    tt.stop();
////    cout<<"Increase BP-Flag update time: "<<tt.GetRuntime()<<" s."<<endl;
////    tt.start();
////    BPFlagUpdateDecrease(overlayBatchDec, true);
//////    BPFlagUpdateDecrease(overlayBatchDec, false);
//            VPLBPFlagConstruct(true);
//            tt.stop();
//            cout<<"Decrease BP-Flag update time: "<<tt.GetRuntime()<<" s."<<endl;

            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::MixRepair_PartiIndexVPLPost, this, true, boost::ref(partiBatchDec), boost::ref(partiBatchInc), boost::ref(tPost) ));
            thread.add_thread(new boost::thread(&Graph::MixRefreshCrossLabelsVPL, this, true, boost::ref(tCross), threadnum-max(affectedParti.size(),affectedPartiInc.size())));
            thread.join_all();
            if(tPost<tCross){
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tPost;
                stageDurations[PH2H_Post]=stageDurations[PH2H_Post]+tCross-tPost;
            }else{
                cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tCross;
//            stageDurations[PH2H_Post]=0;
            }
        }
    }





//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

//function for VPL's throughput test of decrease update
void Graph::VPLBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                if(newW<Neighbor[a][i].second){
                    Neighbor[a][i].second=newW;
                }
                else{
                    cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                if(newW<Neighbor[b][i].second){
                    Neighbor[b][i].second=newW;
                }
                else{
                    cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }

                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(PartiTag[a].second && PartiTag[b].second){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(!PartiTag[a].second && PartiTag[b].second){//b is boundary
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid1) == partiBatch.end()){
                    partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid1].emplace_back(wBatch[k]);
            }else if(PartiTag[a].second && !PartiTag[b].second){//a is boundary
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid2) == partiBatch.end()){
                    partiBatch.insert({pid2,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid2].emplace_back(wBatch[k]);
            }
            else{
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    if(partiBatch.find(pid1) == partiBatch.end()){
                        partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                    }
                    partiBatch[pid1].emplace_back(wBatch[k]);
                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }

            }
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();
    for(int i=0;i<partiNum;++i){
        ProBeginVertexSetParti[i].clear(); vertexIDChLParti[i].clear();
        ProBeginVertexSetPartiExtend[i].clear();
    }
    for(int i=0;i<Tree.size();++i){
        Tree[i].DisRe.clear();
        Tree[i].DisRePost.clear();
    }
    // Step 1: Update shortcuts
    // Approach 1: partitioned version
    // Update the shortcuts of non-boundary vertices
    DecreasePartiBatchUpdateCheckVPL(partiBatch, overlayBatch, true);//multi-thread
//    DecreasePartiBatchUpdateCheckVPL(partiBatch, overlayBatch, false);//single-thread
//    int id1=216209, id2=217962;
//    cout<<"after partition update "<<id1<<"("<<NodeOrder[id1]<<","<<PartiTags[id1].first<<") "<<id2<<"("<<NodeOrder[id2]<<","<<PartiTags[id2].first<<") "<<QueryPostMHL(id1,id2)<<"("<<Dijkstra(id1,id2,Neighbor)<<")"<<endl;
//    for(auto it=Tree[rank[id1]].vert.begin();it!=Tree[rank[id1]].vert.end();++it){
//        if(it->first==id2){
//            cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//            break;
//        }
//    }
    cout<<"overlayBatch size: "<<overlayBatch.size()<<endl;

    // Overlay shortcut update
    DecreaseOverlayBatchVPL(overlayBatch,Tree,rank,heightMax,false);

    // Approach 2: non-partitioned version
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,true);//direct bottom-up, with label construction
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,false);//direct bottom-up, without label construction
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]+=tt2.GetRuntime();
//    cout<<"algoQuery: PCH-No"<<endl;


    vector<int> vParti;
    for(int pid=0;pid<partiNum;++pid){
        if(!ProBeginVertexSetPartiExtend[pid].empty()){
            vParti.push_back(pid);
        }
    }


    if(threadnum<=affectedParti.size()){
        cout<<"No enough thread for parallelization between post-boundary and cross-boundary."<<endl;
        //unparalleled version
        if(algoUpdate>=PH2H_Post){
            double tOverlay=0,tPost=0,tCross=0;
            tt2.start();
            // Step 2: Update overlay index
            DecreaseOverlayBatchLabelVPL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            tt2.stop();
            tOverlay=tt2.GetRuntime();
            cout<<"Overlay index update: "<<tOverlay<<endl;
            tt2.start();
            // Step 3: Update partition index
            Repair_PartiIndexVPLPost(true, false, partiBatch, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, false, partiBatch);//post
            tt2.stop();
            stageDurations[PCH_No]=tt2.GetRuntime();
            if(algoUpdate==PH2H_Cross){
                tt2.start();
                RefreshCrossLabelsVPL(true, false, tCross, threadnum);//extend, parallel
//            RefreshExtensionLabelsPostMHL(false, false);//extend, single-thread
                tt2.stop();
                stageDurations[PH2H_Post]=tt2.GetRuntime();
            }
        }
    }else{//if thread number is sufficient
        cout<<"Post-boundary and cross-boundary could be parallelized!!!"<<endl;
        //parallel version
        double tOverlay=0,tPost=0,tCross=0;
        if(algoUpdate==PH2H_Post) {
            tt2.start();
            DecreaseOverlayBatchLabelVPL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            overlayUpdateT+=tOverlay;
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            tt2.start();
            Repair_PartiIndexVPLPost(true, false, partiBatch, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, false, partiBatch);//post
            tt2.stop();
            stageDurations[PCH_No] = stageDurations[PCH_No]+ tt2.GetRuntime() + tOverlay;
        }else if(algoUpdate==PH2H_Cross){
            tt2.start();
            DecreaseOverlayBatchLabelPostMHL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            overlayUpdateT+=tOverlay;
            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::Repair_PartiIndexVPLPost, this, true, false, boost::ref(partiBatch), boost::ref(tPost) ));
            thread.add_thread(new boost::thread(&Graph::RefreshCrossLabelsVPL, this, true, false, boost::ref(tCross), threadnum-affectedParti.size()));
            thread.join_all();
            if(tPost<tCross){
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tPost;
                stageDurations[PH2H_Post]=stageDurations[PH2H_Post]+tCross-tPost;
            }else{
                cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tCross;
//            stageDurations[PH2H_Post]=0;
            }
        }
    }




    tt.stop();
//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

//function for VPL's throughput test of increase update, optimized version
void Graph::VPLBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    tAncestor=0, tBoundary=0;
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(PartiTag[a].second && PartiTag[b].second){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(!PartiTag[a].second && PartiTag[b].second){
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid1) == partiBatch.end()){
                    partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid1].emplace_back(wBatch[k]);
            }else if(PartiTag[a].second && !PartiTag[b].second){
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid2) == partiBatch.end()){
                    partiBatch.insert({pid2,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid2].emplace_back(wBatch[k]);
            }else
            {
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    if(partiBatch.find(pid1) == partiBatch.end()){
                        partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                    }
                    partiBatch[pid1].emplace_back(wBatch[k]);
                }else{
                    cout<<"Wrong for this edge update. "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<")"<<endl; exit(1);
                }

            }
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    ProBeginVertexSetOverlayInc.clear();
    for(int i=0;i<partiNum;++i){
        ProBeginVertexSetPartiInc[i].clear();
        ProBeginVertexSetPartiExtendInc[i].clear();
    }
    for(int i=0;i<Tree.size();++i){
        Tree[i].DisRe.clear();
        Tree[i].DisRePost.clear();
    }
    // Step 1: Update shortcuts
    // Approach 1: partitioned version
    IncreasePartiBatchUpdateCheckVPL(partiBatch, overlayBatch, true);//multi-thread, bottom-up shortcut update
//    IncreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, false);//single-thread
    cout<<"overlayBatch size: "<<overlayBatch.size()<<endl;
    IncreaseOverlayBatchVPL(overlayBatch,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);

    // Approach 2: non-partitioned version
//    IncreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);//direct bottom-up, with label construction
//    IncreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);//direct bottom-up, without label construction
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]+=tt2.GetRuntime();
//    cout<<"algoQuery: PCH-No"<<endl;

    if(threadnum<=affectedPartiInc.size()){
        // repair the partition index, unparalleled version
        if(algoUpdate>=PH2H_Post){
            tt2.start();
            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlayInc, VidtoTNid);
            tt2.stop();
            double toverlay=tt2.GetRuntime();
            cout<<"Overlay label update time: "<<toverlay<<" s."<<endl;
//        for(int i=0;i<Tree[rank[67975]].disPost.size();++i){
//            if(Tree[rank[67975]].vAncestor[i]==67432){
//                cout<<"B Cnt: "<<Tree[rank[67975]].disPost[i]<<" "<<Tree[rank[67975]].cntPost[i]<<endl;
//            }
//        }
            double tPost=0,tCross=0;
            tt2.start();
            Repair_PartiIndexVPLPost(true, true, partiBatch, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, true, partiBatch);//post
            tt2.stop();
            cout<<"Post-boundary index update time: "<<tt2.GetRuntime()<<" s; Boundary array update: "<<tBoundary<<" s; Ancestor array update: "<<tAncestor<<" s."<<endl;
            stageDurations[PCH_No]=stageDurations[PCH_No]+tt2.GetRuntime()+toverlay;
            if(algoUpdate==PH2H_Cross){
                tt2.start();
                RefreshCrossLabelsVPL(true, true, tCross, threadnum);//extend, single-thread
//            RefreshExtensionLabelsPostMHL(false, true);//extend, single-thread
                algoQuery=PH2H_Cross;
                tt2.stop();
                stageDurations[PH2H_Post]+=tt2.GetRuntime();
            }
        }
    }else{
        double tOverlay=0,tPost=0,tCross=0;
        if(algoUpdate==PH2H_Post) {
            tt2.start();
            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlayInc, VidtoTNid);
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            overlayUpdateT+=tOverlay;
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            tt2.start();
            Repair_PartiIndexVPLPost(true, true, partiBatch, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, true, partiBatch);//post
            tt2.stop();
            stageDurations[PCH_No] = stageDurations[PCH_No]+tt2.GetRuntime() + tOverlay;
        }else if(algoUpdate==PH2H_Cross){
            tt2.start();
            IncreaseOverlayBatchLabelVPL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            overlayUpdateT+=tOverlay;
            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::Repair_PartiIndexVPLPost, this, true, true, boost::ref(partiBatch), boost::ref(tPost) ));
            thread.add_thread(new boost::thread(&Graph::RefreshCrossLabelsVPL, this, true, true, boost::ref(tCross), threadnum-affectedPartiInc.size()));
            thread.join_all();
            if(tPost<tCross){
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tPost;
                stageDurations[PH2H_Post]=stageDurations[PH2H_Post]+tCross-tPost;
            }else{
                cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tCross;
//            stageDurations[PH2H_Post]=0;
            }
        }
    }



    tt.stop();
//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;

//    CorrectnessCheck(100);
}

//function for throughput test of decrease update
void Graph::PostMHLBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    tt2.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                if(newW<Neighbor[a][i].second){
                    Neighbor[a][i].second=newW;
                }
                else{
                    cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                if(newW<Neighbor[b][i].second){
                    Neighbor[b][i].second=newW;
                }
                else{
                    cout<<"Wrong decrease update. "<<a<<" "<<b<<" "<< Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }

                break;
            }
        }

        int pid1=PartiTags[a].first, pid2=PartiTags[b].first;
        if(pid1==-1 && pid2==-1){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(pid1!=-1 && pid2==-1){
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid1) == partiBatch.end()){
                    partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid1].emplace_back(wBatch[k]);
            }else if(pid1==-1 && pid2!=-1){
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid2) == partiBatch.end()){
                    partiBatch.insert({pid2,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid2].emplace_back(wBatch[k]);
            }else
            {
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    if(partiBatch.find(pid1) == partiBatch.end()){
                        partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                    }
                    partiBatch[pid1].emplace_back(wBatch[k]);
                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }

            }
        }
    }
    tt2.stop();
    uStageDurations[Dijk]+=tt2.GetRuntime();
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();
    for(int i=0;i<partiNum;++i){
        ProBeginVertexSetParti[i].clear(); vertexIDChLParti[i].clear();
        ProBeginVertexSetPartiExtend[i].clear();
    }
    for(int i=0;i<Tree.size();++i){
        Tree[i].DisRe.clear();
        Tree[i].DisRePost.clear();
    }
    // Step 1: Update shortcuts
    // Approach 1: partitioned version
    DecreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, true);//multi-thread
//    DecreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, false);//single-thread
//    int id1=216209, id2=217962;
//    cout<<"after partition update "<<id1<<"("<<NodeOrder[id1]<<","<<PartiTags[id1].first<<") "<<id2<<"("<<NodeOrder[id2]<<","<<PartiTags[id2].first<<") "<<QueryPostMHL(id1,id2)<<"("<<Dijkstra(id1,id2,Neighbor)<<")"<<endl;
//    for(auto it=Tree[rank[id1]].vert.begin();it!=Tree[rank[id1]].vert.end();++it){
//        if(it->first==id2){
//            cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//            break;
//        }
//    }
    cout<<"overlayBatch size: "<<overlayBatch.size()<<endl;
    DecreaseOverlayBatchPostMHL(overlayBatch,Tree,rank,heightMax,false);
//    cout<<"after overlay update "<<id1<<"("<<NodeOrder[id1]<<","<<PartiTags[id1].first<<") "<<id2<<"("<<NodeOrder[id2]<<","<<PartiTags[id2].first<<") "<<QueryPostMHL(id1,id2)<<"("<<Dijkstra(id1,id2,Neighbor)<<")"<<endl;
//    for(auto it=Tree[rank[id1]].vert.begin();it!=Tree[rank[id1]].vert.end();++it){
//        if(it->first==id2){
//            cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//            break;
//        }
//    }
    // Approach 2: non-partitioned version
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,true);//direct bottom-up, with label construction
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,false);//direct bottom-up, without label construction
    algoQuery=PCH_No;
    tt2.stop();
    uStageDurations[PCH_No]+=tt2.GetRuntime();
    stageDurations[Dijk]+=tt2.GetRuntime();
//    cout<<"algoQuery: PCH-No"<<endl;

    vector<int> vParti;
    for(int pid=0;pid<partiNum;++pid){
        if(!ProBeginVertexSetPartiExtend[pid].empty()){
            vParti.push_back(pid);
        }
    }


    if(threadnum<=affectedParti.size()){
        cout<<"No enough thread for parallelization between post-boundary and cross-boundary."<<endl;
        //unparalleled version
        if(algoUpdate>=PH2H_Post){
            double tOverlay=0,tPost=0,tCross=0;
            tt2.start();
            // Step 2: Update overlay index
            DecreaseOverlayBatchLabelPostMHL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            tt2.stop();
            tOverlay=tt2.GetRuntime();
            cout<<"Overlay index update: "<<tOverlay<<endl;
            uStageDurations[PH2H_No] += tOverlay;
            tt2.start();
            // Step 3: Update partition index
            Repair_PartiIndexPostMHLPost(true, DecUpdate, partiBatch, tPost);//post
            tt2.stop();
            stageDurations[PCH_No]=stageDurations[PCH_No]+tt2.GetRuntime()+tOverlay;
            uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tt2.GetRuntime();
            cout<<"Post-boundary update time: "<<tPost<< " s"<<endl;
            if(algoUpdate==PH2H_Cross){
                tt2.start();
                RefreshExtensionLabelsPostMHL(true, false, tCross, threadnum);//extend, single-thread
//            RefreshExtensionLabelsPostMHL(false, false);//extend, single-thread
                tt2.stop();
                stageDurations[PH2H_Post]+=tt2.GetRuntime();
                uStageDurations[PH2H_Cross] += tt2.GetRuntime();
                cout<<"Cross-boundary update time: "<<tCross<< " s"<<endl;
            }
        }
    }
    else{//if thread number is sufficient
        cout<<"Post-boundary and cross-boundary could be parallelized!!!"<<endl;
        //parallel version
        double tOverlay=0,tPost=0,tCross=0;
        if(algoUpdate==PH2H_Post) {
            tt2.start();
            DecreaseOverlayBatchLabelPostMHL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            uStageDurations[PH2H_No] += tOverlay;
            overlayUpdateT+=tOverlay;
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            tt2.start();
            Repair_PartiIndexPostMHLPost(true, DecUpdate, partiBatch, tPost);//post
            tt2.stop();
            stageDurations[PCH_No] = stageDurations[PCH_No]+ tt2.GetRuntime() + tOverlay;
            uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tt2.GetRuntime();
        }else if(algoUpdate==PH2H_Cross){
            tt2.start();
            DecreaseOverlayBatchLabelPostMHL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            uStageDurations[PH2H_No] += tOverlay;
            overlayUpdateT+=tOverlay;
            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::Repair_PartiIndexPostMHLPost, this, true, DecUpdate, boost::ref(partiBatch), boost::ref(tPost) ));
            thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelsPostMHL, this, true, false, boost::ref(tCross), threadnum-2));//threadnum-affectedParti.size()
            thread.join_all();
            if(tPost<tCross){
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tPost;
                stageDurations[PH2H_Post]=stageDurations[PH2H_Post]+tCross-tPost;
                uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tPost;
                uStageDurations[PH2H_Cross]=uStageDurations[PH2H_Cross]+tCross-tPost;
            }else{
                cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tCross;
                uStageDurations[PH2H_Cross] = uStageDurations[PH2H_Cross]+tCross;
//            stageDurations[PH2H_Post]=0;
            }
            cout<<"Post and Cross-boundary update time: "<<tPost<<" "<<tCross<< " s"<<endl;
        }
    }




    tt.stop();
//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

//function for throughput test of increase update, optimized version
void Graph::PostMHLBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    tAncestor=0, tBoundary=0;
    Timer tt;
    Timer tt2;
    tt.start();
    tt2.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTags[a].first, pid2=PartiTags[b].first;
        if(pid1==-1 && pid2==-1){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(pid1!=-1 && pid2==-1){
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid1) == partiBatch.end()){
                    partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid1].emplace_back(wBatch[k]);
            }else if(pid1==-1 && pid2!=-1){
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid2) == partiBatch.end()){
                    partiBatch.insert({pid2,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid2].emplace_back(wBatch[k]);
            }else
            {
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    if(partiBatch.find(pid1) == partiBatch.end()){
                        partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                    }
                    partiBatch[pid1].emplace_back(wBatch[k]);
                }else{
                    cout<<"Wrong for this edge update. "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<")"<<endl; exit(1);
                }

            }
        }
    }
    tt2.stop();
    algoQuery=Dijk;
    uStageDurations[Dijk]+=tt2.GetRuntime();
    tt2.start();
    vUpdated.assign(node_num, false);
    ProBeginVertexSetOverlay.clear();
    for(int i=0;i<partiNum;++i){
        ProBeginVertexSetParti[i].clear();
        ProBeginVertexSetPartiExtend[i].clear();
    }
    for(int i=0;i<Tree.size();++i){
        Tree[i].DisRe.clear();
        Tree[i].DisRePost.clear();
    }
    // Step 1: Update shortcuts
    // Approach 1: partitioned version
    IncreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, true);//multi-thread, bottom-up shortcut update
//    IncreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, false);//single-thread
    cout<<"overlayBatch size: "<<overlayBatch.size()<<endl;
    IncreaseOverlayBatchPostMHL(overlayBatch,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);

    // Approach 2: non-partitioned version
//    IncreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);//direct bottom-up, with label construction
//    IncreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);//direct bottom-up, without label construction
    algoQuery=PCH_No;
    tt2.stop();
    uStageDurations[PCH_No]+=tt2.GetRuntime();
    stageDurations[Dijk]+=tt2.GetRuntime();
//    cout<<"algoQuery: PCH-No"<<endl;

    if(threadnum<=affectedParti.size()){
        // repair the partition index, unparalleled version
        if(algoUpdate>=PH2H_Post){
            tt2.start();
            IncreaseOverlayBatchLabelPostMHL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
            tt2.stop();
            double toverlay=tt2.GetRuntime();
            cout<<"Overlay label update time: "<<toverlay<<" s."<<endl;
            uStageDurations[PH2H_No] += toverlay;
            double tPost=0,tCross=0;
            tt2.start();
            Repair_PartiIndexPostMHLPost(true, IncUpdate, partiBatch, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, true, partiBatch);//post
            tt2.stop();
            cout<<"Post-boundary index update time: "<<tt2.GetRuntime()<<" s; Boundary array update: "<<tBoundary<<" s; Ancestor array update: "<<tAncestor<<" s."<<endl;
            stageDurations[PCH_No]=stageDurations[PCH_No]+tt2.GetRuntime()+toverlay;
            uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tt2.GetRuntime();
            if(algoUpdate==PH2H_Cross){
                tt2.start();
                RefreshExtensionLabelsPostMHL(true, true, tCross, threadnum);//extend, single-thread
//            RefreshExtensionLabelsPostMHL(false, true);//extend, single-thread
                algoQuery=PH2H_Cross;
                tt2.stop();
                uStageDurations[PH2H_Cross] += tt2.GetRuntime();
                stageDurations[PH2H_Post]+=tt2.GetRuntime();
            }
        }
    }else{
        double tOverlay=0,tPost=0,tCross=0;
        if(algoUpdate==PH2H_Post) {
            tt2.start();
            IncreaseOverlayBatchLabelPostMHL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            overlayUpdateT+=tOverlay;
            uStageDurations[PH2H_No] += tOverlay;
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            tt2.start();
            Repair_PartiIndexPostMHLPost(true, IncUpdate, partiBatch, tPost);//post
            tt2.stop();
            stageDurations[PCH_No] = stageDurations[PCH_No]+tt2.GetRuntime() + tOverlay;
            uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tt2.GetRuntime();
        }else if(algoUpdate==PH2H_Cross){
            tt2.start();
            IncreaseOverlayBatchLabelPostMHL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
            tt2.stop();
            tOverlay = tt2.GetRuntime();
            uStageDurations[PH2H_No] += tOverlay;
            cout << "Overlay label update time: " << tOverlay << " s." << endl;
            overlayUpdateT+=tOverlay;
            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::Repair_PartiIndexPostMHLPost, this, true, IncUpdate, boost::ref(partiBatch), boost::ref(tPost) ));
            thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelsPostMHL, this, true, true, boost::ref(tCross), threadnum-2));//threadnum-affectedParti.size()
            thread.join_all();
            if(tPost<tCross){
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tPost;
                stageDurations[PH2H_Post]=stageDurations[PH2H_Post]+tCross-tPost;
                uStageDurations[PH2H_Post] = uStageDurations[PH2H_Post]+tPost;
                uStageDurations[PH2H_Cross]=uStageDurations[PH2H_Cross]+tCross-tPost;
            }else{
                cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
                stageDurations[PCH_No] = stageDurations[PCH_No]+tOverlay+tCross;
                uStageDurations[PH2H_Cross] = uStageDurations[PH2H_Cross]+tCross;
//            stageDurations[PH2H_Post]=0;
            }
        }
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;

//    CorrectnessCheck(100);
}

//function for throughput test of increase update
void Graph::PMHLBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());

    if(!partiBatch.empty()){
        if(threadnum==1){
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
                IncreasePartiBatchUpdateCheckCH(pid, it->second, overlayBatch, false, updatedSCs[pid]);
            }
        }else{
            if(partiBatch.size()>threadnum){
                cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
                vector<vector<int>> processID;
                processID.assign(threadnum, vector<int>());
                vector<int> vertices;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    vertices.emplace_back(pid);
                }
                ThreadDistribute(vertices, processID);
                boost::thread_group thread;
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheckCHV, this, processID[j], boost::ref(partiBatch), boost::ref(overlayBatch), false, boost::ref(updatedSCs) ));
                }
                thread.join_all();
            }else{
                boost::thread_group thread;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false, boost::ref(updatedSCs[pid]) ));
                }
                thread.join_all();
            }
        }

    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;

        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }
//            cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
            if(newdis>olddis){//if '=', not problem; if '<', problem
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//                sm->notify();
            }
            else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
//    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<< updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);
    }

    cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]+=tt2.GetRuntime();
    tt2.start();


    if(threadnum==1){
        IncreaseOverlayBatchLabel(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
        if(!partiBatch.empty()){
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
                IncreasePartiBatchLabel(Trees[pid], ranks[pid], heightMaxs[pid], ProBeginVertexSetParti[pid], VidtoTNidP);//vector<Node> &Tree, vector<int> &rank, int heightMax, int& checknum, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid
            }
        }
    }else{
        if(partiBatch.size()+1>threadnum){
            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::IncreaseOverlayBatchLabel, this, boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(ProBeginVertexSetOverlay), boost::ref(VidtoTNid) ));
            vector<vector<int>> processID;
            processID.assign(threadnum-1, vector<int>());
            vector<int> vertices;
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchLabelV, this,  processID[j], boost::ref(Trees), boost::ref(ranks), boost::ref(heightMaxs), boost::ref(ProBeginVertexSetParti), boost::ref(VidtoTNidP) ));
            }
            thread.join_all();
        }else{
            boost::thread_group thread;
            thread.add_thread(new boost::thread(&Graph::IncreaseOverlayBatchLabel, this, boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(ProBeginVertexSetOverlay), boost::ref(VidtoTNid) ));
            if(!partiBatch.empty()){
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchLabel, this, boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], boost::ref(ProBeginVertexSetParti[pid]), boost::ref(VidtoTNidP) ));//vector<Node> &Tree, vector<int> &rank, int heightMax, int& checknum, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid
                }
            }
            thread.join_all();
        }
    }


    algoQuery=PH2H_No;
    tt2.stop();
    stageDurations[PCH_No]+=tt2.GetRuntime();
    // repair the partition index
    if(algoUpdate>=PH2H_Post){
//        Trees=TreesNo;
        tt2.start();
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);
        algoQuery=PH2H_Post;
        tt2.stop();
        stageDurations[PH2H_No]+=tt2.GetRuntime();
        if(algoUpdate==PH2H_Cross){
            tt2.start();
//            RefreshExtensionLabels(partiBatch);
            if(threadnum==1){
                RefreshExtensionLabelsNoAllPair(partiBatch, false);
            }else{
                RefreshExtensionLabelsNoAllPair(partiBatch, true);
            }

            algoQuery=PH2H_Cross;
            tt2.stop();
            stageDurations[PH2H_Post]+=tt2.GetRuntime();
        }
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_No]+stageDurations[PH2H_Post];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;

//    CorrectnessCheck(100);
}

//function for throughput test of increase update, optimized version
void Graph::PMHLBatchUpdateIncOpt(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
//            cout<<"Overlay update: "<<a<<"("<<PartiTag[a].first<<","<<PartiTag[a].second<<") "<<b<<"("<<PartiTag[b].first<<","<<PartiTag[b].second<<") "<<oldW<<" "<<NeighborsOverlay[a][b]<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
//            cout<<"Parti update: "<<pid1<<". "<<a<<"("<<PartiTag[a].first<<","<<PartiTag[a].second<<") "<<" "<<b<<"("<<PartiTag[b].first<<","<<PartiTag[b].second<<") "<<" "<<oldW<<" "<<newW<<endl;
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
//    cout<<"Before parti CH, size of overlayBatch: "<<overlayBatch.size()<<endl;
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), true, boost::ref(updatedSCs[pid]) ));
//            IncreasePartiBatchUpdateCheckCH(pid, it->second, overlayBatch, true );
        }
        thread.join_all();

    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;

        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }
//        cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
            if(newdis>olddis){//if '=', not problem; if '<', problem
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
//                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<<updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);
    }

    cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;

//    cout<<"After parti CH, size of overlayBatch: "<<overlayBatch.size()<<endl;
    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    algoQuery=PCH_No;
    if(algoUpdate>PCH_No) {
        IncreaseOverlayBatchLabel(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
    }

    // repair the partition index
    if(algoUpdate>=PH2H_Post){
        Repair_PartiIndex(true, true, partiBatch);//post
//        Repair_PartiIndexForOpt(false, true, partiBatch);//post

        algoQuery=PH2H_Post;
        if(algoUpdate==PH2H_Cross){
//            RefreshExtensionLabels(partiBatch);
            RefreshExtensionLabelsNoAllPair(partiBatch,true);
            algoQuery=PH2H_Cross;
        }
    }

    tt.stop();
    runT1+=tt.GetRuntime();
    cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;

//    CorrectnessCheck(100);
}

//function for throughput test of decrease update
void Graph::DecBatchThroughput(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    vUpdated.assign(node_num, false);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();
    }

    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,true);

    algoQuery=PH2H_No;
    // repair the partition index
    if(algoUpdate>=PH2H_Post){
        Repair_PartiIndex(true, false, partiBatch);//post
        algoQuery=PH2H_Post;
        if(algoUpdate==PH2H_Cross){
            RefreshExtensionLabels(partiBatch);//extend
            algoQuery=PH2H_Cross;
        }
    }

    tt.stop();
    runT1+=tt.GetRuntime();
    cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
}
//function for throughput test of increase update
void Graph::IncBatchThroughput(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
//            cout<<"Overlay update. "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<oldW<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
//            cout<<"Parti update. "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<oldW<<" "<<newW<<endl;
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    vUpdated.assign(node_num, false);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();

    }

    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);
    algoQuery=PH2H_No;

    // repair the partition index
    if(algoUpdate>=PH2H_Post){
//        Trees=TreesNo;
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);
        algoQuery=PH2H_Post;
        if(algoUpdate==PH2H_Cross){
            RefreshExtensionLabels(partiBatch);
            algoQuery=PH2H_Cross;
        }
    }

    tt.stop();
    runT1+=tt.GetRuntime();
    cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
}


//function for throughput test of decrease update
void Graph::DecBatchThroughputNP(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    Timer tt2;
    tt.start();

    map<int,int> checkedDis;
//    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the fresh distance and avoid search in the adjacent list

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }


//    vector<set<int>> SCre; //SCre.clear(); not truly define the shortcut priority
    vector<set<OrderCompp>> SCre; //SCre.clear(); truly define the shortcut priority
//    SCre.assign(node_num,set<int>());//{vertexID, set<int>}
    SCre.assign(node_num,set<OrderCompp>());
    set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    int a,b,oldW,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second; oldW=wBatch[k].second.first;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
                if(newW<Neighbor[a][i].second){
                    Neighbor[a][i].second=newW;
                }else{
                    cout<<"Invalid update. "<<a<<" "<<b<<" "<<Neighbor[a][i].second<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
                if(newW<Neighbor[b][i].second){
                    Neighbor[b][i].second=newW;
                }
                else{
                    cout<<"Invalid update. "<<a<<" "<<b<<" "<<Neighbor[b][i].second<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }
    }

    algoQuery=Dijk;

    tt2.start();
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second; oldW=wBatch[k].second.first;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }


        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompp(lid));
//                    OCdis[make_pair(lid,hid)]=newW;
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }
    }

    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw;
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(Tree[rank[ProID]].dis[cidH]>Cw){
                Tree[rank[ProID]].dis[cidH]=Cw;
                Tree[rank[ProID]].FN[cidH]=true;
                Tree[rank[ProID]].cnt[cidH]=1;//new
                ProIDdisCha=true;
                Tree[rank[ProID]].DisRe.insert(Cid);
            }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                Tree[rank[ProID]].FN[cidH]=true;
                Tree[rank[ProID]].cnt[cidH]+=1;//new
            }

            int hid,hidHeight,lid,lidHeight,wsum;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                if(Hnei.find(hid)!=Hnei.end()){
                    wsum=Cw+Hnei[hid];
                    if(wsum<Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.first=wsum;
                        Tree[rank[Cid]].vert[j].second.second=1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
//                        OCdis[make_pair(Cid,hid)]=wsum;
                    }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid]].vert[k].second.first>wsum){
                            Tree[rank[lid]].vert[k].second.first=wsum;
                            Tree[rank[lid]].vert[k].second.second=1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
//                            OCdis[make_pair(lid,Cid)]=wsum;
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            if(SCre[ProID].find(lid)==SCre[ProID].end()){//if not found, to avoid repeat count
//                            if(OCdis.find(make_pair(ProID,lid))==OCdis.end()){//if not found
                                Tree[rank[lid]].vert[k].second.second+=1;
                            }
//                            Tree[rank[lid]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChL.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQuery(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }
    }

    algoQuery=CH;
    tt2.stop();
    stageDurations[Dijk]+=tt2.GetRuntime();
    tt2.start();
    //cout<<"Finish bottom-up refresh"<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5H2H(rank[ProBeginVertexID], linee, vertexIDChL,checkedDis);
    }
    //return checkedDis.size();
    algoQuery=H2H;
    tt2.stop();
    stageDurations[CH]+=tt2.GetRuntime();
    tt.stop();
//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[CH];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[CH]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

void Graph::EachNodeProBDis5H2H(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
    bool ProIDdisCha=false;
    vector<int> cntNew(line.size(),0);
    vector<bool> flagUpdate(line.size(),false);
    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];//use cntNew to redress the cnt value since the edge decrease may lead to more ways for dis (i.e., increase the cnt)
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                Tree[child].cnt[i]=1;//new
                                ProIDdisCha=true;
                                flagUpdate[i]=true;
                                cntNew[i]=1;
                            }
                            else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                                cntNew[i]++;
                                if(flagUpdate[i]) {
                                    Tree[child].cnt[i]+=1;//new
                                }
                                else if(cntNew[i]>Tree[child].cnt[i]){
                                    Tree[child].cnt[i]=cntNew[i];
                                }
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }

                }
            }
        }
    }
    else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
                        }
                        else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                        Tree[child].FN[i]=false;
                        Tree[child].cnt[i]=1;//new
                        ProIDdisCha=true;
                        flagUpdate[i]=true;
                        cntNew[i]=1;
                    }
                    else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                        cntNew[i]++;
                        if(flagUpdate[i]) {
                            Tree[child].cnt[i]+=1;//new
                        }
                        else if(cntNew[i]>Tree[child].cnt[i]){
                            Tree[child].cnt[i]=cntNew[i];
                        }
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5H2H(Tree[child].ch[i], line, vertexIDChL,checkedDis);
    }
    line.pop_back();

}

//function for throughput test of decrease update, original one
/*void Graph::IncBatchThroughputNP(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    Timer tt2;
    tt.start();

    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]), original weight of the affected shortcut, maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    //NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear(); the affected shortcut pair
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompp> OC; OC.clear();//the lower-order vertex of the affected shortcut, vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW<newW){
            bool ifFind=false;
            for(int i=0;i<Neighbor[a].size();i++){
                if(Neighbor[a][i].first==b){
                    Neighbor[a][i].second=newW;
                    ifFind=true;
                    break;
                }
            }
            assert(ifFind);
            ifFind=false;
            for(int i=0;i<Neighbor[b].size();i++){
                if(Neighbor[b][i].first==a){
                    Neighbor[b][i].second=newW;
                    ifFind=true;
                    break;
                }
            }
            assert(ifFind);

        }else{
            cout<<"Wrong update."<<oldW<<" "<<newW<<endl; exit(1);
        }
    }

    algoQuery=Dijk;
    tt2.start();
    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW<newW){

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }

            for(int i=0;i<Tree[rank[lid]].vert.size();i++){
                if(Tree[rank[lid]].vert[i].first==hid){
                    if(Tree[rank[lid]].vert[i].second.first==oldW){
                        Tree[rank[lid]].vert[i].second.second-=1;
                        if(Tree[rank[lid]].vert[i].second.second<1){//the shortcut needs update, should be increased
                            OCdis[make_pair(lid,hid)]=oldW;//original weight of the affected shortcut
                            SCre[lid].insert(hid);//the affected shortcut pair
                            OC.insert(OrderCompp(lid));//the lower-order vertex of the affected shortcut
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID; vector<int> line;
    /// Shortcut update
    while(!OC.empty()){
        ProID=(*OC.begin()).x;//from the lowest-order vertex
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        // get the ancestors of ProID, each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[pachid]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.emplace_back(Vert[j].first,Vert[j].second.first);
                }
            }
            //check the affected shortcuts
            int hid,lid;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;
                if(Hnei.find(hid)!=Hnei.end()){
                    if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second-=1;
                        if(Tree[rank[Cid]].vert[j].second.second<1){
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompp(Cid));
                            OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                            Tree[rank[lid]].vert[k].second.second-=1;
                            if(Tree[rank[lid]].vert[k].second.second<1){
                                SCre[lid].insert(Cid);
                                OC.insert(OrderCompp(lid));
                                OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }

//get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;

            for(int i=0;i<Neighbor[ProID].size();i++){
                if(Neighbor[ProID][i].first==Cid){
                    newCw=Neighbor[ProID][i].second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(!Wnodes.empty()){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first;
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

//            if(newCw<=Cw){
//                cout<<"Invalid update 1. "<<ProID<<" "<<Cid<<" "<<Cw<<"("<<newCw<<") ";
//                for(int i=0;i<Tree[rank[ProID]].vert.size();i++) {
//                    if (Tree[rank[ProID]].vert[i].first == Cid) {
//                        cout<<Tree[rank[ProID]].vert[i].second.first<<"("<<newCw<<") "<<Tree[rank[ProID]].vert[i].second.second<<"("<<countwt<<")";
//                    }
//                }
//                cout<<endl;
//            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                if(Tree[rank[ProID]].vert[i].first==Cid){
                    Tree[rank[ProID]].vert[i].second.first=newCw;
                    Tree[rank[ProID]].vert[i].second.second=countwt;
                    break;
                }
            }

            if(newCw>Cw){
                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProID]].FN[cidH]){//if the distance label is from shortcut, then the label may be affected.
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProID]].FN[cidH]=false;//? may still be the source
                    Tree[rank[ProID]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }
                }
            }



        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQuery(rnew,r)!=rnew){//if they are in different branches
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;//identify the roots
        }

    }
    algoQuery=CH;
    tt2.stop();
    stageDurations[Dijk]=tt2.GetRuntime();

    if(algoUpdate==H2H){
        tt2.start();
        int ProBeginVertexID;
//    cout<<"Root number: "<<ProBeginVertexSet.size()<<endl;
        for(int i=0;i<ProBeginVertexSet.size();i++){//for each root
            ProBeginVertexID=ProBeginVertexSet[i];
//        cout<<i<<" "<<ProBeginVertexID<<" "<<Tree[rank[ProBeginVertexID]].height<<endl;
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1H2H(rank[ProBeginVertexID], linee,checknum);
        }
        //return checknum;
        algoQuery=H2H;
        tt2.stop();
        stageDurations[CH]=tt2.GetRuntime();
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
    runT1=runT1+stageDurations[Dijk]+stageDurations[CH];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[CH]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}*/

//function for throughput test of decrease update, new
void Graph::IncBatchThroughputNP(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    Timer tt2;
    tt.start();

    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]), original weight of the affected shortcut, maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    //NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<OrderCompp>> SCre; //SCre.clear(); the affected shortcut pair
    SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
    set<OrderCompp> OC; OC.clear();//the lower-order vertex of the affected shortcut, vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW<newW){
            bool ifFind=false;
            for(int i=0;i<Neighbor[a].size();i++){
                if(Neighbor[a][i].first==b){
                    Neighbor[a][i].second=newW;
                    ifFind=true;
                    break;
                }
            }
            assert(ifFind);
            ifFind=false;
            for(int i=0;i<Neighbor[b].size();i++){
                if(Neighbor[b][i].first==a){
                    Neighbor[b][i].second=newW;
                    ifFind=true;
                    break;
                }
            }
            assert(ifFind);

        }else{
            cout<<"Wrong update."<<oldW<<" "<<newW<<endl; exit(1);
        }
    }

    algoQuery=Dijk;
    tt2.start();
    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW<newW){

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }

            for(int i=0;i<Tree[rank[lid]].vert.size();i++){
                if(Tree[rank[lid]].vert[i].first==hid){
                    if(Tree[rank[lid]].vert[i].second.first==oldW){
                        Tree[rank[lid]].vert[i].second.second-=1;
                        if(Tree[rank[lid]].vert[i].second.second<1){//the shortcut needs update, should be increased
                            OCdis[make_pair(lid,hid)]=oldW;//original weight of the affected shortcut
                            SCre[lid].insert(hid);//the affected shortcut pair
                            OC.insert(OrderCompp(lid));//the lower-order vertex of the affected shortcut
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID; vector<int> line;
    /// Shortcut update
    while(!OC.empty()){
        ProID=(*OC.begin()).x;//from the lowest-order vertex
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        // get the ancestors of ProID, each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[pachid]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.emplace_back(Vert[j].first,Vert[j].second.first);
                }
            }
            //check the affected shortcuts
            int hid,lid;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;
                if(Hnei.find(hid)!=Hnei.end()){
                    if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second-=1;
                        if(Tree[rank[Cid]].vert[j].second.second<1){
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompp(Cid));
                            OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                            if(SCre[ProID].find(lid)==SCre[ProID].end()) {//if not found, to avoid repeat count
                                Tree[rank[lid]].vert[k].second.second -= 1;
                                if (Tree[rank[lid]].vert[k].second.second < 1) {
                                    SCre[lid].insert(Cid);
                                    OC.insert(OrderCompp(lid));
                                    OCdis[make_pair(lid, Cid)] = Cw + Lnei[j].second;
                                }
                            }
                        }
                        break;
                    }
                }
            }

//get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;

            for(int i=0;i<Neighbor[ProID].size();i++){
                if(Neighbor[ProID][i].first==Cid){
                    newCw=Neighbor[ProID][i].second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(!Wnodes.empty()){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first;
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

            if(newCw<=Cw){
                cout<<"Invalid update 1. "<<ProID<<" "<<Cid<<" "<<Cw<<"("<<newCw<<") ";
                for(int i=0;i<Tree[rank[ProID]].vert.size();i++) {
                    if (Tree[rank[ProID]].vert[i].first == Cid) {
                        cout<<Tree[rank[ProID]].vert[i].second.first<<"("<<newCw<<") "<<Tree[rank[ProID]].vert[i].second.second<<"("<<countwt<<")";
                    }
                }
                cout<<endl;
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                if(Tree[rank[ProID]].vert[i].first==Cid){
                    Tree[rank[ProID]].vert[i].second.first=newCw;
                    Tree[rank[ProID]].vert[i].second.second=countwt;
                    break;
                }
            }

            if(newCw>Cw){
                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProID]].FN[cidH]){//if the distance label is from shortcut, then the label may be affected.
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProID]].FN[cidH]=false;//? may still be the source
                    Tree[rank[ProID]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }
                }
            }



        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQuery(rnew,r)!=rnew){//if they are in different branches
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;//identify the roots
        }

    }
    algoQuery=CH;
    tt2.stop();
    stageDurations[Dijk]+=tt2.GetRuntime();

    if(algoUpdate==H2H){
        tt2.start();
        int ProBeginVertexID;
//    cout<<"Root number: "<<ProBeginVertexSet.size()<<endl;
        for(int i=0;i<ProBeginVertexSet.size();i++){//for each root
            ProBeginVertexID=ProBeginVertexSet[i];
//        cout<<i<<" "<<ProBeginVertexID<<" "<<Tree[rank[ProBeginVertexID]].height<<endl;
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1H2H(rank[ProBeginVertexID], linee,checknum);
        }
        //return checknum;
        algoQuery=H2H;
        tt2.stop();
        stageDurations[CH]+=tt2.GetRuntime();
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
//    runT1=runT1+stageDurations[Dijk]+stageDurations[CH];
//    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[CH]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

void Graph::eachNodeProcessIncrease1H2H(int children, vector<int>& line, int& changelabel){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
    for(int i=0;i<Tree[children].dis.size();i++){
        if(Tree[children].cnt[i]<=0){//if the distance label to i-th ancestor should be maintained, cnt may less than actual value
//        if(true){
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int PID;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){//for the tree node that contains childID as vert element
                PID=VidtoTNid[childID][k];
                if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){//if label is from shortcut
                    Tree[PID].cnt[i]-=1;
                }
            }

            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                PID=VidtoTNid[line[i]][k];
//                if(Tree[PID].height>Tree[children].height){///modified for correctness, PID may not be the descendant of children
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){///modified for correctness, PID may not be the descendant of children
                    if(PID>Tree.size()){
                        cout<<"PID error! "<<PID<<" "<<Tree.size()<<endl; exit(1);
                    }
                    if(childH>Tree[PID].dis.size()){
                        cout<<"childH error! "<<childH<<" "<<Tree[PID].dis.size()<<": "<<children<<"("<<Tree[children].height<<") "<<PID<<"("<<Tree[PID].height<<")"<<endl; exit(1);
                    }
                    if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){///
                        Tree[PID].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[b]].height-1;
                if(bH<i){
                    if(Dvb+Tree[rank[line[i]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[line[i]]].dis[bH];
                        count=1;
                    }else if(Dvb+Tree[rank[line[i]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{
                    if(Dvb+Tree[rank[b]].dis[i]<dis){
                        dis=Dvb+Tree[rank[b]].dis[i];
                        count=1;
                    }else if(Dvb+Tree[rank[b]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) {
                Tree[children].FN[i]=true;
            }
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1H2H(Tree[children].ch[i],line,changelabel);
    }
    line.pop_back();
}

//function for efficiency test
unsigned long long Graph::EffiCheckThroughput(vector<pair<int,int>>& ODpair, Timer& tRecord, int batchInterval, unsigned long long& throughputNum) {
    bool ifDebug=false;
//    ifDebug=true;
    int s, t;
    double runT = 0, runT0=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
    int d1, d2;
    Timer tt;
    unsigned long long runtimes = 0, runtimes0=0, runtimes1=0, runtimes2=0, runtimes3=0, runtimes4=0, runtimes5=0;
    vector<int> results(ODpair.size(), -1);
    int i = 0;
    double tNow = 0;
    if(algoChoice==1){//CH+H2H
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == CH) {//CH
                tt.start();
                d2 = QueryCHWP(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == H2H) {//H2H
                tt.start();
                d2 = QueryH2H(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT2 += tt.GetRuntime();
                ++runtimes2;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery==Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
                }
            }

        }

        cout<<"CH+H2H. Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (CH): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (H2H): Duration: "<<runT2<<" s; Throughput number: "<<runtimes2<<" ; query time: "<<1000 * runT2 / runtimes2 << " ms."<<endl;
    }
    else if(algoChoice==2){//PH2H
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == PH2H_No) {//PH2H-No
                tt.start();
                d2 = Query(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == PH2H_Post) {//PH2H-Post
                tt.start();
                d2 = Query(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT2 += tt.GetRuntime();
                ++runtimes2;
            }
            else if (algoQuery == PH2H_Cross) {//PH2H-Extend
                tt.start();
                d2 = Query(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT3 += tt.GetRuntime();
                ++runtimes3;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery==Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" ("<<PartiTag[s].first<<", "<<PartiTag[s].second<<") "<<t<<" ("<<PartiTag[t].first<<", "<<PartiTag[t].second<<") "<<d2<<" "<<d1<<endl;
                    exit(1);
                }
            }
        }

        cout<<"PH2H. Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (No-boundary): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (Post-boundary): Duration: "<<runT2<<" s; Throughput number: "<<runtimes2<<" ; query time: "<<1000 * runT2 / runtimes2 << " ms."<<endl;
        cout<<"Stage 4 (Extension): Duration: "<<runT3<<" s; Throughput number: "<<runtimes3<<" ; query time: "<<1000 * runT3 / runtimes3 << " ms."<<endl;
    }
    else if(algoChoice==3){//PCH+PH2H
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == PCH_No) {//PCH-No
                tt.start();
                d2 = QueryPMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == PH2H_No) {//PH2H-No
                tt.start();
                d2 = QueryPMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT2 += tt.GetRuntime();
                ++runtimes2;
            }
            else if (algoQuery == PH2H_Post) {//PH2H-Post
                tt.start();
                d2 = QueryPMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT4 += tt.GetRuntime();
                ++runtimes4;
            }
            else if (algoQuery == PH2H_Cross) {//PH2H-Extend
                tt.start();
                d2 = QueryPMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT5 += tt.GetRuntime();
                ++runtimes5;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery== Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" ("<<PartiTag[s].first<<", "<<PartiTag[s].second<<") "<<t<<" ("<<PartiTag[t].first<<", "<<PartiTag[t].second<<") "<<d2<<" "<<d1<<endl;
                    exit(1);
                }
            }
        }

        cout<<"PCH+PH2H. Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (PCH-No): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-No): Duration: "<<runT2<<" s; Throughput number: "<<runtimes2<<" ; query time: "<<1000 * runT2 / runtimes2 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Post): Duration: "<<runT4<<" s; Throughput number: "<<runtimes4<<" ; query time: "<<1000 * runT4 / runtimes4 << " ms."<<endl;
        cout<<"Stage 5 (PH2H-Extend): Duration: "<<runT5<<" s; Throughput number: "<<runtimes5<<" ; query time: "<<1000 * runT5 / runtimes5 << " ms."<<endl;
    }
    else if(algoChoice==4){//PCH+PH2H with optimization
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == PCH_No) {//PCH-No
                tt.start();
                d2 = QueryPMHLOpt(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == PH2H_Post) {//PCH-Post
                tt.start();
                d2 = QueryPMHLOpt(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT3 += tt.GetRuntime();
                ++runtimes3;
            }
            else if (algoQuery == PH2H_Post) {//PH2H-Post
                tt.start();
                d2 = QueryPMHLOpt(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT4 += tt.GetRuntime();
                ++runtimes4;
            }
            else if (algoQuery == PH2H_Cross) {//PH2H-Extend
                tt.start();
                d2 = QueryPMHLOpt(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT5 += tt.GetRuntime();
                ++runtimes5;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery== Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" ("<<PartiTag[s].first<<", "<<PartiTag[s].second<<") "<<t<<" ("<<PartiTag[t].first<<", "<<PartiTag[t].second<<") "<<d2<<" "<<d1<<endl;
                    exit(1);
                }
            }
        }

        cout<<"PCH+PH2H (Opt). Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (PCH-No): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-Post): Duration: "<<runT4<<" s; Throughput number: "<<runtimes4<<" ; query time: "<<1000 * runT4 / runtimes4 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Extend): Duration: "<<runT5<<" s; Throughput number: "<<runtimes5<<" ; query time: "<<1000 * runT5 / runtimes5 << " ms."<<endl;
    }
    else if(algoChoice==5){//PostMHL
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == PCH_No) {//PCH-No
                tt.start();
                d2 = QueryPostMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == PH2H_Post) {//PH2H-Post
                tt.start();
                d2 = QueryPostMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT4 += tt.GetRuntime();
                ++runtimes4;
            }
            else if (algoQuery == PH2H_Cross) {//PH2H-Extend
                tt.start();
                d2 = QueryPostMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT5 += tt.GetRuntime();
                ++runtimes5;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery== Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" ("<<PartiTags[s].first<<") "<<t<<" ("<<PartiTags[t].first<<") "<<d2<<" "<<d1<<endl;
                    QueryPostMHLDebug(s,t);
                    exit(1);
                }
            }
        }

        cout<<"PostMHL. Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (PCH-No): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-Post): Duration: "<<runT4<<" s; Throughput number: "<<runtimes4<<" ; query time: "<<1000 * runT4 / runtimes4 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Extend): Duration: "<<runT5<<" s; Throughput number: "<<runtimes5<<" ; query time: "<<1000 * runT5 / runtimes5 << " ms."<<endl;
    }

    throughputNum+=runtimes;
    return runtimes;
}

void Graph::IndexMaintenance(int updateType, int updateSize, bool ifBatch, int batchSize) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = sourcePath+"/"+dataset + ".update";
    bool ifDebug=false;
    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand(0);
    int updateBatch=1;
    updateBatch=max(updateBatch,updateSize/batchSize);
    cout<<"Update batch: "<<updateBatch<<" ; Batch size: "<<batchSize<<endl;
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);

    Timer tt;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 0:{
            break;
        }
        case 1:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;
            Graph g2=*this;
            if(ifBatch){//for batch update
                if(updateBatch*batchSize>updateData.size()){
                    updateBatch=floor(updateData.size()/batchSize);
                }
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<updateBatch;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        update_i=u*batchSize+i;
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;
                        oldW = updateData[update_i].second;
                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

//                    g2.QueryPartiPartiExtDebug(3321,212184);

                    tt.start();
                    g2.DecreaseBatch(wBatch);
//                    DecreaseBatch(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease update Time: "<<runT1/(updateBatch*batchSize)<<" s.\n"<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<updateBatch;u++){
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*0.5;
                    if(newW < 1) {
                        cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        exit(1);
                    }
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    }
                    tt.start();
                    g2.DecreaseSingle(ID1,ID2,oldW,newW);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }

                }

                cout<<"Average Decrease update Time: "<<runT1/updateBatch<<" s.\n"<<endl;
            }

//            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            ifIncrease=true;
            if(ifBatch){//for batch update
                if(updateBatch*batchSize>updateData.size()){
                    updateBatch=floor(updateData.size()/batchSize);
                }
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=1;u<updateBatch;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        update_i=u*batchSize+i;
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;
                        oldW = updateData[update_i].second;
                        newW=oldW*1.5;
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
//                        ++update_i;
                    }

//                    QueryCoreDebug(111009,104090);
//                    QueryPartiPartiExtDebug(107830,104726);

                    tt.start();
                    IncreaseBatch(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }
                cout<<"Average Increase update Time: "<<runT2/(updateBatch*batchSize)<<" s.\n"<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<updateBatch;u++){
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*2;
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<"("<<PartiTag[ID1].first<<") "<<ID2<<"("<<PartiTag[ID2].first<<") "<<oldW<<" "<<newW<<endl;
                    }

                    tt.start();
                    IncreaseSingle(ID1,ID2,oldW,newW);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }

                cout<<"Average Increase update Time: "<<runT2/updateBatch<<" s.\n"<<endl;
            }

            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

void Graph::DecreaseBatch(vector<pair<pair<int, int>, pair<int, int>>> &wBatch) {
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }
        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();
    }

    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,true);

    // repair the partition index
    if(algoUpdate>=PH2H_Post){
        Repair_PartiIndex(true, false, partiBatch);//post
        if(algoUpdate==PH2H_Cross){
            RefreshExtensionLabels(partiBatch);//extend
        }
    }
}
//partition update of PH2H
void Graph::DecreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay){
    //partition batch decrease update
    vector<pair<pair<int,int>,int>> updatedSC;

    DecreasePartiBatch(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, true);

    vector<int> Bid=BoundVertex[pid];
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(int i=0;i<Bid.size();i++){
        bid1=Bid[i];
        for(int j=i+1;j<Bid.size();j++){
            bid2=Bid[j];
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                continue;//exit(1);
            }

            newdis=QueryH2HPartition(bid1,bid2,pid);
            if(newdis<olddis){
//                    cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
                sm->wait();
                weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
                sm->notify();
            }
        }
    }
}
//partition update PCH
void Graph::DecreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC){
    //partition batch decrease update
//    vector<pair<pair<int,int>,int>> updatedSC;


    if(ifOpt){
        DecreasePartiBatchForOpt(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, false, false);
    }else{
        DecreasePartiBatch(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, false);
    }


//    vector<int> Bid=BoundVertex[pid];
//    //check the boundary edge within partition
//    int bid1,bid2,olddis,newdis;
//    for(auto it=updatedSC.begin();it!=updatedSC.end();++it){
//        bid1=it->first.first, bid2=it->first.second; newdis=it->second;
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }
//
//        if(newdis<olddis){
////            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
////                NeighborsOverlay[bid1][bid2]=newdis;
////                NeighborsOverlay[bid2][bid1]=newdis;
//            sm->wait();
//            weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
//        }
//    }

}

//partition update PCH
void Graph::DecreasePartiBatchUpdateCheckCHV(vector<int> p, map<int,vector<pair<pair<int,int>,pair<int,int>>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<vector<pair<pair<int,int>,int>>>& updatedSC){
    //partition batch decrease update
//    vector<pair<pair<int,int>,int>> updatedSC;

    for(int i=0;i<p.size();++i){
        int pid=p[i];
        if(ifOpt){
            DecreasePartiBatchForOpt(pid, wBatch[pid], NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC[pid], false, false);
        }else{
            DecreasePartiBatch(pid, wBatch[pid], NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC[pid], false);
        }
    }
}

//shortcut update for PostMHL
void Graph::DecreasePartiBatchUpdateCheckPostMHL(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int,int>,pair<int,int>>>& overlayBatch, bool ifParallel){
    //partition batch decrease update
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());

    if(ifParallel){
        if(!partiBatch.empty()){
            if(partiBatch.size()>threadnum){
                cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
                vector<vector<int>> processID;
                processID.assign(threadnum, vector<int>());
                vector<int> vertices;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    vertices.emplace_back(pid);
                }
                ThreadDistribute(vertices, processID);
                boost::thread_group thread;
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchPostMHLShortcutV, this, processID[j], boost::ref(partiBatch), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCs) ));
                }
                thread.join_all();

            }
            else{
//            cout<<"Update Partition number: "<<partiBatch.size()<<endl;
                boost::thread_group thread;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchPostMHLShortcut, this, pid, boost::ref(it->second), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCs[pid]) ));
                }
                thread.join_all();
            }

        }
    }else{
        cout<<"Single thread computation."<<endl;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
//            cout<<"Shortcut update of partition "<<pid<<endl;
            DecreasePartiBatchPostMHLShortcut(pid, it->second, Tree, rank, heightMax, updatedSCs[pid]);
        }
    }


//    cout<<"updatedSC size: "<<updatedSC.size()<<endl;
    map<pair<int,int>, int> updatedSCTrue;
    int updateSCSize=0;
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second, newdis=it->second;
            updateSCSize++;
            if(updatedSCTrue.find(make_pair(bid1,bid2))==updatedSCTrue.end()){//if not found
                updatedSCTrue.insert({make_pair(bid1,bid2),newdis});
            }else {//if found
//            cout<<bid1<<" "<<bid2<<": "<<updatedSCTrue[make_pair(bid1,bid2)]<<" "<<newdis<<endl;
                if(updatedSCTrue[make_pair(bid1,bid2)]>newdis){//if found
                    updatedSCTrue[make_pair(bid1,bid2)]=newdis;
                }
            }
        }
    }

//    cout<<"updatedSCTrue size: "<<updatedSCTrue.size()<<" "<<updateSCSize<< endl;

    int sVertex;
    for(auto it=updatedSCTrue.begin();it!=updatedSCTrue.end();++it){
        bid1=it->first.first, bid2=it->first.second, newdis=it->second;
//        for(int i=0;i<SCconNodesMT[bid1][bid2].size();i++){
//            sVertex=SCconNodesMT[bid1][bid2][i].first;
//            if(NeighborsOverlay[bid1][bid2] > Tree[rank[sVertex]].)
//            SCconNodesMT[ID1][ID2].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
//        }

        olddis=-1;
        for(auto it2=Tree[rank[bid1]].vert.begin();it2!=Tree[rank[bid1]].vert.end();++it2){
            if(it2->first==bid2){
                olddis=it2->second.first;
                break;
            }
        }
        if(olddis==-1){
            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
        }
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }

        if(newdis<=olddis){//cnt may change
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//            sm->wait();
            overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
        }else if(newdis>olddis){
            cout<<"Something wrong happens. "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl; exit(1);
        }
    }

}

// Decrease shortcut update for VPL
void Graph::DecreasePartiBatchUpdateCheckVPL(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int,int>,pair<int,int>>>& overlayBatch, bool ifParallel){
    //partition batch decrease update
    if(!partiBatch.empty()){
        vector<vector<pair<pair<int,int>,int>>> updatedSCs;
        updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());

        if(ifParallel & partiBatch.size()>1){
            if(partiBatch.size()>threadnum){
                cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
                vector<vector<int>> processID;
                processID.assign(threadnum, vector<int>());
                vector<int> vertices;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    vertices.emplace_back(pid);
                }
                ThreadDistribute(vertices, processID);
                boost::thread_group thread;
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchVPLShortcutV, this, processID[j], boost::ref(partiBatch), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCs) ));
                }
                thread.join_all();

            }
            else{
//            cout<<"Update Partition number: "<<partiBatch.size()<<endl;
                boost::thread_group thread;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchVPLShortcut, this, pid, boost::ref(it->second), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCs[pid]) ));
                }
                thread.join_all();
            }
        }
        else{
//            cout<<"Single thread computation."<<endl;
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
//            cout<<"Shortcut update of partition "<<pid<<endl;
                DecreasePartiBatchVPLShortcut(pid, it->second, Tree, rank, heightMax, updatedSCs[pid]);
            }
        }


//    cout<<"updatedSC size: "<<updatedSC.size()<<endl;
        map<pair<int,int>, int> updatedSCTrue;
        int updateSCSize=0;
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
            int pid=it1->first;
            for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
                bid1=it->first.first, bid2=it->first.second, newdis=it->second;
                updateSCSize++;
                if(updatedSCTrue.find(make_pair(bid1,bid2))==updatedSCTrue.end()){//if not found
                    updatedSCTrue.insert({make_pair(bid1,bid2),newdis});
                }else {//if found
//            cout<<bid1<<" "<<bid2<<": "<<updatedSCTrue[make_pair(bid1,bid2)]<<" "<<newdis<<endl;
                    if(updatedSCTrue[make_pair(bid1,bid2)]>newdis){//if found
                        updatedSCTrue[make_pair(bid1,bid2)]=newdis;
                    }
                }
            }
        }

//    cout<<"updatedSCTrue size: "<<updatedSCTrue.size()<<" "<<updateSCSize<< endl;

        int sVertex;
        for(auto it=updatedSCTrue.begin();it!=updatedSCTrue.end();++it){
            bid1=it->first.first, bid2=it->first.second, newdis=it->second;
//        for(int i=0;i<SCconNodesMT[bid1][bid2].size();i++){
//            sVertex=SCconNodesMT[bid1][bid2][i].first;
//            if(NeighborsOverlay[bid1][bid2] > Tree[rank[sVertex]].)
//            SCconNodesMT[ID1][ID2].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
//        }

            olddis=-1;
            for(auto it2=Tree[rank[bid1]].vert.begin();it2!=Tree[rank[bid1]].vert.end();++it2){
                if(it2->first==bid2){
                    olddis=it2->second.first;
                    break;
                }
            }
            if(olddis==-1){
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }

            if(newdis<=olddis){//cnt may change
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//            sm->wait();
                overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
            }else if(newdis>olddis){
                cout<<"dec Something wrong happens. "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl; exit(1);//not for mix
            }
        }
    }
}


void Graph::IncreaseBatch(vector<pair<pair<int, int>, pair<int, int>>> &wBatch) {
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }
        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();

    }

    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);

    // repair the partition index
    if(algoUpdate>=PH2H_Post){
//        Trees=TreesNo;
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);
        if(algoUpdate==PH2H_Cross){
            RefreshExtensionLabels(partiBatch);
        }
    }
}
//for PH2H
void Graph::IncreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay){
    //partition batch Increase update
    vector<pair<pair<int,int>,int>> updatedSC;

    IncreasePartiBatch(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,true);

    vector<int> Bid=BoundVertex[pid];
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(int i=0;i<Bid.size();i++){
        bid1=Bid[i];
        for(int j=i+1;j<Bid.size();j++){
            bid2=Bid[j];
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                continue;//exit(1);
            }

            newdis=QueryH2HPartition(bid1,bid2,pid);
            if(newdis>olddis){//if '=', not problem; if '<', problem
                sm->wait();
                weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
}
//for PCH
void Graph::IncreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC){
    //partition batch Increase update
//    vector<pair<pair<int,int>,int>> updatedSC;

    if(ifOpt){
        IncreasePartiBatchForOpt(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,false);
    }else{
        IncreasePartiBatch(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,false);
    }

//    cout<<pid<<". Size of updatedSC: "<<updatedSC.size()<<endl;
//    vector<int> Bid=BoundVertex[pid];
//    //check the boundary edge within partition
//    int bid1,bid2,olddis,newdis;
//
//    for(auto it=updatedSC.begin();it!=updatedSC.end();++it){
//        bid1=it->first.first, bid2=it->first.second; newdis=it->second;
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }
////        cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//        if(newdis>olddis){//if '=', not problem; if '<', problem
//            sm->wait();
//            weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//            sm->notify();
//        } else if(newdis<olddis){
//            cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
//            exit(1);
//        }
//    }

}

void Graph::IncreasePartiBatchUpdateCheckCHV(vector<int> p, map<int,vector<pair<pair<int,int>,pair<int,int>>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& overlayBatch, bool ifOpt, vector<vector<pair<pair<int,int>,int>>>& updatedSCs){
    for(int i=0;i<p.size();++i){
        int pid=p[i];
        IncreasePartiBatchUpdateCheckCH(pid, wBatch[pid], overlayBatch, false, updatedSCs[pid] );
    }
}

//bottom-up shortcut update for post-boundary partition index
void Graph::IncreasePartiBatchUpdateCheckPostMHL(map<int, vector<pair<pair<int, int>, pair<int, int>>>> &partiBatch, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatch, bool ifParallel) {
    //partition batch Increase update
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;//the overlay shortcuts that may be affected
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
    vector<int> PropagateOverlay;

    if(ifParallel){
        if(!partiBatch.empty()){
            if(partiBatch.size()>threadnum){
                cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
                vector<vector<int>> processID;
                processID.assign(threadnum, vector<int>());
                vector<int> vertices;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    vertices.emplace_back(pid);
                }
                ThreadDistribute(vertices, processID);
                boost::thread_group thread;
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchPostMHLShortcutV, this, processID[j], boost::ref(partiBatch), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCs), boost::ref(PropagateOverlay) ));
                }
                thread.join_all();

            }
            else{
//            cout<<"Update Partition number: "<<partiBatch.size()<<endl;
                boost::thread_group thread;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchPostMHLShortcut, this, pid, boost::ref(it->second), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCs[pid]), boost::ref(PropagateOverlay) ));
                }
                thread.join_all();
            }

        }
    }else{
        cout<<"Single thread computation."<<endl;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
//            cout<<"Shortcut update of partition "<<pid<<endl;
            IncreasePartiBatchPostMHLShortcut(pid, it->second, Tree, rank, heightMax, updatedSCs[pid], PropagateOverlay);
        }
    }

//    cout<<"updatedSC size: "<<updatedSC.size()<<endl;
    map<pair<int,int>, pair<int,int>> updatedSCTrue;
    int updateSCSize=0;
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis,cnt, wsum;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second, wsum=it->second;
            updateSCSize++;
//            cout<<bid1<<" "<<bid2<<": "<<newdis<<endl;
            if(updatedSCTrue.find(make_pair(bid1,bid2))==updatedSCTrue.end()){//if not found
                updatedSCTrue.insert({make_pair(bid1,bid2),make_pair(wsum,1)});
            }else {//if found
//            cout<<bid1<<" "<<bid2<<": "<<updatedSCTrue[make_pair(bid1,bid2)]<<" "<<newdis<<endl;
                if(updatedSCTrue[make_pair(bid1,bid2)].first!=wsum){
                    cout<<"1. Inconsistent wsum! "<<bid1<<" "<<bid2<<" "<<updatedSCTrue[make_pair(bid1,bid2)].first<<" "<<wsum<<endl; exit(1);
                }
                updatedSCTrue[make_pair(bid1,bid2)].second+=1;
            }
        }
    }

//    cout<<"updatedSCTrue size: "<<updatedSCTrue.size()<<" "<<updateSCSize<< endl;
//    cout<<"Initial overlayBatch size: "<<weightOverlay.size()<<endl;
    int sVertex;
    bool ifFound;
    for(auto it=updatedSCTrue.begin();it!=updatedSCTrue.end();++it){
        bid1=it->first.first, bid2=it->first.second, wsum=it->second.first, cnt=it->second.second;
        ifFound=false;
        for(int i=0;i<Tree[rank[bid1]].vert.size();i++){
            if(Tree[rank[bid1]].vert[i].first==bid2){
                ifFound=true;
                if(Tree[rank[bid1]].vert[i].second.second < 1){//need update
                    newdis=INF; int countwt=0;
                    for(auto it2=Neighbor[bid1].begin();it2!=Neighbor[bid1].end();++it2){
                        if(it2->first==bid2){
                            newdis=it2->second;//the weight value in the original graph
                            countwt=1;
                            break;
                        }
                    }
                    int ssw,wtt,wid;
                    vector<pair<int,int>> Wnodes;
                    Wnodes.clear();
                    if(bid1<bid2)
                        Wnodes=SCconNodesMT[bid1][bid2]; //cout<<"wid num "<<Wnodes.size()<<endl;
                    else
                        Wnodes=SCconNodesMT[bid2][bid1];
                    if(Wnodes.size()>0){
                        for(int i=0;i<Wnodes.size();i++){
                            wid=Wnodes[i].first;
                            for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                                if(Tree[rank[wid]].vert[j].first==bid1){
                                    ssw=Tree[rank[wid]].vert[j].second.first;
                                }
                                if(Tree[rank[wid]].vert[j].first==bid2){
                                    wtt=Tree[rank[wid]].vert[j].second.first;
                                }
                            }

                            if(ssw+wtt<newdis){
                                newdis=ssw+wtt;
                                countwt=1;
                            }else if(ssw+wtt==newdis){
                                countwt+=1;
                            }
                        }
                    }

//                    olddis=-1;
//                    for(auto it2=NeighborsOverlayV[bid1].begin();it2!=NeighborsOverlayV[bid1].end();++it2){
//                        if(it2->first==bid2){
//                            olddis=it2->second;
//                            break;
//                        }
//                    }
//                    if(olddis==-1){
//                        cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//                    }
//                    if(olddis!=wsum){
//                        cout<<"wsum is inconsitent with olddis! "<<bid1<<" "<<bid2<<" "<<wsum<<" "<<olddis<<endl; exit(1);//？
//                    }
                    olddis=wsum;
                    if(newdis>olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//            sm->wait();
                        overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
                    }else if(newdis<olddis){
                        cout<<"Something wrong happens. "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl; exit(1);
                    }
                }
                else{
//                    cout<<"No update for this boundary shortcut! "<<bid1<<" "<<bid2<<" "<<Tree[rank[bid1]].vert[i].second.first<<" "<<Tree[rank[bid1]].vert[i].second.second<<" "<<cnt<<endl;
                }
                break;
            }
        }
        if(!ifFound){
            cout<<"Not found "<<bid2<<" in "<<bid1<<" 's vert. "<<endl ;exit(1);
        }
    }
    //get the LCA of PropagateOverlay
//    if(!PropagateOverlay.empty()){
////        cout<<"PropagateOverlay size: "<<PropagateOverlay.size()<<endl;
//        int ProID=PropagateOverlay[0];
//        int rLCA;
//        int VID;
//        for(int i=1;i<PropagateOverlay.size();++i){
//            VID=PropagateOverlay[i];
//            rLCA = LCAQueryOverlay(rank[ProID], rank[VID]);
//            ProID = Tree[rLCA].uniqueVertex;
//        }
//        ProBeginVertexSetOverlay.push_back(ProID);
//    }


}


void Graph::MixVPLBatchShortcutUpdateParti(map<int, vector<pair<pair<int, int>, pair<int, int>>>> &partiBatchDec, map<int, vector<pair<pair<int, int>, pair<int, int>>>> &partiBatchInc,vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchDec, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatchInc, bool ifParallel){
    map<pair<int, int>, pair<int, int>> updatedSCTrueInc;
    map<pair<int,int>, int> updatedSCTrueDec;
    int bid1,bid2,olddis,newdis,wsum,cnt;

    //partition batch mix update
    if(!partiBatchDec.empty() || !partiBatchInc.empty()){
        vector<vector<pair<pair<int,int>,int>>> updatedSCsDec;
        updatedSCsDec.assign(partiNum,vector<pair<pair<int,int>,int>>());
        vector<vector<pair<pair<int,int>,int>>> updatedSCsInc;
        updatedSCsInc.assign(partiNum,vector<pair<pair<int,int>,int>>());
        map<int, pair<vector<pair<pair<int, int>, pair<int, int>>>,vector<pair<pair<int, int>, pair<int, int>>>>> partiBatch;
        for(auto it=partiBatchDec.begin();it!=partiBatchDec.end();++it){
            partiBatch[it->first].first=it->second;
        }
        for(auto it=partiBatchInc.begin();it!=partiBatchInc.end();++it){
            partiBatch[it->first].second=it->second;
        }
//        if(false){
        if(ifParallel & partiBatch.size()>1){
            if(partiBatch.size()>threadnum){
                cout<<"partiBatch is larger than thread number! "<<partiBatchDec.size()<<" "<<threadnum<< endl;
                vector<vector<int>> processID;
                processID.assign(threadnum, vector<int>());
                vector<int> vertices;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    vertices.emplace_back(pid);
                }
                ThreadDistribute(vertices, processID);
                boost::thread_group thread;
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::MixPartiBatchVPLShortcutV, this, boost::ref(processID[j]), boost::ref(partiBatch), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCsDec),  boost::ref(updatedSCsInc) ));
                }
                thread.join_all();

            }
            else{
//            cout<<"Update Partition number: "<<partiBatch.size()<<endl;
                boost::thread_group thread;
                for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                    int pid=it->first;
                    thread.add_thread(new boost::thread(&Graph::MixPartiBatchVPLShortcut, this, pid, boost::ref(it->second.first), boost::ref(it->second.second), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCsDec[pid]),  boost::ref(updatedSCsInc[pid]) ));
                }
                thread.join_all();
            }
        }
        else{
//            cout<<"Single thread computation."<<endl;
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
//            cout<<"Shortcut update of partition "<<pid<<endl;
                MixPartiBatchVPLShortcut(pid, it->second.first, it->second.second, Tree, rank, heightMax, updatedSCsDec[pid], updatedSCsInc[pid]);
            }
        }


//    cout<<"updatedSC size: "<<updatedSC.size()<<endl;

        int updateSCSize=0;
        //check the boundary edge within partition
        for(auto it=overlayBatchDec.begin();it!=overlayBatchDec.end();++it){
            if(NodeOrder[it->first.first]<NodeOrder[it->first.second]){
                bid1=it->first.first, bid2=it->first.second;
            }else{
                bid1=it->first.second, bid2=it->first.first;
            }
            newdis=it->second.second;
            if(updatedSCTrueDec.find(make_pair(bid1,bid2))==updatedSCTrueDec.end()){//if not found
                updatedSCTrueDec.insert({make_pair(bid1,bid2),newdis});
            }else {//if found
                cout<<"Wrong. "<<bid1<<" "<<bid2<<": "<<updatedSCTrueDec[make_pair(bid1,bid2)]<<" "<<newdis<<endl; exit(1);
            }
            updateSCSize++;
        }
        for(auto it1=partiBatchDec.begin();it1!=partiBatchDec.end();++it1){
            int pid=it1->first;
            for(auto it=updatedSCsDec[pid].begin();it!=updatedSCsDec[pid].end();++it){
                bid1=it->first.first, bid2=it->first.second, newdis=it->second;
                updateSCSize++;
                if(updatedSCTrueDec.find(make_pair(bid1,bid2))==updatedSCTrueDec.end()){//if not found
                    updatedSCTrueDec.insert({make_pair(bid1,bid2),newdis});
                }else {//if found
//                    cout<<pid<<" "<<bid1<<" "<<bid2<<": "<<updatedSCTrueDec[make_pair(bid1,bid2)]<<" "<<newdis<<endl;
                    if(updatedSCTrueDec[make_pair(bid1,bid2)]>newdis){//if found
                        updatedSCTrueDec[make_pair(bid1,bid2)]=newdis;
                    }
                }
            }
        }

//        cout<<"updatedSCTrueDec size: "<<updatedSCTrueDec.size()<<" "<<updateSCSize<< endl;
        updateSCSize=0;
        for(auto it=overlayBatchInc.begin();it!=overlayBatchInc.end();++it){
            if(NodeOrder[it->first.first]<NodeOrder[it->first.second]){
                bid1=it->first.first, bid2=it->first.second;
            }else{
                bid1=it->first.second, bid2=it->first.first;
            }

            olddis=it->second.first, newdis=it->second.second;
            if(updatedSCTrueInc.find(make_pair(bid1,bid2))==updatedSCTrueInc.end()){//if not found
                updatedSCTrueInc.insert({make_pair(bid1,bid2),make_pair(olddis,newdis)});
            }else {//if found
                cout<<"Wrong. "<<bid1<<" "<<bid2<<": "<<updatedSCTrueInc[make_pair(bid1,bid2)].second<<" "<<newdis<<endl; exit(1);
            }
            updateSCSize++;
        }
        for (auto it1 = partiBatchInc.begin(); it1 != partiBatchInc.end(); ++it1) {
            int pid = it1->first;
            for (auto it = updatedSCsInc[pid].begin(); it != updatedSCsInc[pid].end(); ++it) {
                bid1 = it->first.first, bid2 = it->first.second, wsum = it->second;
                updateSCSize++;
//            cout<<bid1<<" "<<bid2<<": "<<newdis<<endl;
                if (updatedSCTrueInc.find(make_pair(bid1, bid2)) == updatedSCTrueInc.end()) {//if not found
                    updatedSCTrueInc.insert({make_pair(bid1, bid2), make_pair(wsum, 1)});
                } else {//if found
//            cout<<bid1<<" "<<bid2<<": "<<updatedSCTrue[make_pair(bid1,bid2)]<<" "<<newdis<<endl;
//                    updatedSCTrueInc[make_pair(bid1, bid2)].first = wsum;
                    if (updatedSCTrueInc[make_pair(bid1, bid2)].first != wsum) {
                        cout << "1. Inconsistent wsum! " << bid1<<"("<<PartiTag[bid1].first<<","<<PartiTag[bid1].second << ") " << bid2 << "("<<PartiTag[bid2].first<<","<<PartiTag[bid2].second << ") " <<  updatedSCTrueInc[make_pair(bid1, bid2)].first << " " << wsum << endl;
                        if (updatedSCTrueInc[make_pair(bid1, bid2)].first > wsum) {
                            updatedSCTrueInc[make_pair(bid1, bid2)].first = wsum;
                        }
//                        exit(1);
                    }
                    updatedSCTrueInc[make_pair(bid1, bid2)].second += 1;
                }
            }
        }
//        cout<<"updatedSCTrueInc size: "<<updatedSCTrueInc.size()<<" "<<updateSCSize<< endl;

    }

    if(!updatedSCTrueInc.empty()){
        int sVertex;
        for (auto it = updatedSCTrueInc.begin(); it != updatedSCTrueInc.end(); ++it) {
            bid1 = it->first.first, bid2 = it->first.second, wsum = it->second.first, cnt = it->second.second;
//            if(bid1==613236 && bid2==623280){
//                cout<<"Find inc. "<<bid1<<" "<<bid2<<" "<<wsum<<" "<<cnt<<endl;
//            }
//            if(updatedSCTrueDec.find(it->first)!=updatedSCTrueDec.end()){//if found in updatedSCTrueInc
//                cout<<"Find dec in Inc. "<<bid1<<" "<<bid2<<" "<<wsum<<" "<<cnt<<" "<<updatedSCTrueDec[it->first]<<endl;
//                exit(1);
//                continue;
//            }
            bool ifFound = false;
            for (int i = 0; i < Tree[rank[bid1]].vert.size(); i++) {
                if (Tree[rank[bid1]].vert[i].first == bid2) {
                    ifFound = true;
                    if (Tree[rank[bid1]].vert[i].second.second < 1) {//need update
                        newdis = INF;
                        int countwt = 0;
                        for (auto it2 = Neighbor[bid1].begin(); it2 != Neighbor[bid1].end(); ++it2) {
                            if (it2->first == bid2) {
                                newdis = it2->second;//the weight value in the original graph
                                countwt = 1;
                                break;
                            }
                        }
                        int ssw, wtt, wid;
                        vector<pair<int, int>> Wnodes;
                        Wnodes.clear();
                        if (bid1 < bid2)
                            Wnodes = SCconNodesMT[bid1][bid2]; //cout<<"wid num "<<Wnodes.size()<<endl;
                        else
                            Wnodes = SCconNodesMT[bid2][bid1];
                        if (Wnodes.size() > 0) {
                            for (int i = 0; i < Wnodes.size(); i++) {
                                wid = Wnodes[i].first;
                                for (int j = 0; j < Tree[rank[wid]].vert.size(); j++) {
                                    if (Tree[rank[wid]].vert[j].first == bid1) {
                                        ssw = Tree[rank[wid]].vert[j].second.first;
                                    }
                                    if (Tree[rank[wid]].vert[j].first == bid2) {
                                        wtt = Tree[rank[wid]].vert[j].second.first;
                                    }
                                }

                                if (ssw + wtt < newdis) {
                                    newdis = ssw + wtt;
                                    countwt = 1;
                                } else if (ssw + wtt == newdis) {
                                    countwt += 1;
                                }
                            }
                        }

                        olddis = wsum;
                        if (newdis > olddis) {
                            overlayBatchInc.emplace_back(make_pair(bid1, bid2), make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
                        } else if (newdis < olddis) {
                            overlayBatchDec.emplace_back(make_pair(bid1,bid2), make_pair(olddis,newdis));
//                            cout << "find decrease shortcut update in inc. " << bid1 << " " << bid2 << " " << wsum<<"("<<cnt<<") "<< olddis << " " << newdis <<"("<<Dijkstra(bid1,bid2,Neighbor)<<")"<< endl;//not for mix
//                            exit(1);
                        }
                    } else {
//                    cout<<"No update for this boundary shortcut! "<<bid1<<" "<<bid2<<" "<<Tree[rank[bid1]].vert[i].second.first<<" "<<Tree[rank[bid1]].vert[i].second.second<<" "<<cnt<<endl;
                    }
                    break;
                }
            }
            if (!ifFound) {
                cout << "Not found " << bid2 << " in " << bid1 << " 's vert. " << endl;
                exit(1);
            }
        }
    }

    if(!updatedSCTrueDec.empty()){
        int sVertex;
        for(auto it=updatedSCTrueDec.begin();it!=updatedSCTrueDec.end();++it){
            bid1=it->first.first, bid2=it->first.second, newdis=it->second;
//            if(bid1==613236 && bid2==623280){
//                cout<<"Find dec. "<<bid1<<" "<<bid2<<" "<<newdis<<endl;
//            }
            olddis=-1;
            for(auto it2=Tree[rank[bid1]].vert.begin();it2!=Tree[rank[bid1]].vert.end();++it2){
                if(it2->first==bid2){
                    olddis=it2->second.first;
                    break;
                }
            }
            if(olddis==-1){
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }

            if(newdis<=olddis){//cnt may change
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
                overlayBatchDec.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
            }else if(newdis>olddis){
                cout<<"dec Something wrong happens. "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                exit(1);//not for mix
            }

        }
    }
}

//bottom-up shortcut update for post-boundary partition index
void Graph::IncreasePartiBatchUpdateCheckVPL(map<int, vector<pair<pair<int, int>, pair<int, int>>>> &partiBatch, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatch, bool ifParallel) {
    //partition batch Increase update
    if(!partiBatch.empty()) {
        vector<vector<pair<pair<int, int>, int>>> updatedSCs;//the overlay shortcuts that may be affected
        updatedSCs.assign(partiNum, vector<pair<pair<int, int>, int>>());
        vector<int> PropagateOverlay;

        if (ifParallel & partiBatch.size()>1) {
            if (partiBatch.size() > threadnum) {
                cout << "partiBatch is larger than thread number! " << partiBatch.size() << " " << threadnum << endl;
                vector<vector<int>> processID;
                processID.assign(threadnum, vector<int>());
                vector<int> vertices;
                for (auto it = partiBatch.begin(); it != partiBatch.end(); ++it) {
                    int pid = it->first;
                    vertices.emplace_back(pid);
                }
                ThreadDistribute(vertices, processID);
                boost::thread_group thread;
                for (auto j = 0; j < processID.size(); ++j) {
                    thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchVPLShortcutV, this, processID[j],
                                                        boost::ref(partiBatch), boost::ref(Tree), boost::ref(rank),
                                                        heightMax, boost::ref(updatedSCs)));
                }
                thread.join_all();

            } else {
//            cout<<"Update Partition number: "<<partiBatch.size()<<endl;
                boost::thread_group thread;
                for (auto it = partiBatch.begin(); it != partiBatch.end(); ++it) {
                    int pid = it->first;
                    thread.add_thread(
                            new boost::thread(&Graph::IncreasePartiBatchVPLShortcut, this, pid, boost::ref(it->second),
                                              boost::ref(Tree), boost::ref(rank), heightMax,
                                              boost::ref(updatedSCs[pid])));
                }
                thread.join_all();
            }


        }
        else {
//            cout << "Single thread computation." << endl;
            for (auto it = partiBatch.begin(); it != partiBatch.end(); ++it) {
                int pid = it->first;
//            cout<<"Shortcut update of partition "<<pid<<endl;
                IncreasePartiBatchVPLShortcut(pid, it->second, Tree, rank, heightMax, updatedSCs[pid]);
            }
        }

//    cout<<"updatedSC size: "<<updatedSC.size()<<endl;
        map<pair<int, int>, pair<int, int>> updatedSCTrue;
        int updateSCSize = 0;
        //check the boundary edge within partition
        int bid1, bid2, olddis, newdis, cnt, wsum;
        for (auto it1 = partiBatch.begin(); it1 != partiBatch.end(); ++it1) {
            int pid = it1->first;
            for (auto it = updatedSCs[pid].begin(); it != updatedSCs[pid].end(); ++it) {
                bid1 = it->first.first, bid2 = it->first.second, wsum = it->second;
                updateSCSize++;
//            cout<<bid1<<" "<<bid2<<": "<<newdis<<endl;
                if (updatedSCTrue.find(make_pair(bid1, bid2)) == updatedSCTrue.end()) {//if not found
                    updatedSCTrue.insert({make_pair(bid1, bid2), make_pair(wsum, 1)});
                } else {//if found
//            cout<<bid1<<" "<<bid2<<": "<<updatedSCTrue[make_pair(bid1,bid2)]<<" "<<newdis<<endl;
                    if (updatedSCTrue[make_pair(bid1, bid2)].first != wsum) {
                        cout << "1. Inconsistent wsum! " << bid1 << " " << bid2 << " "
                             << updatedSCTrue[make_pair(bid1, bid2)].first << " " << wsum << endl;
                        exit(1);
                    }
                    updatedSCTrue[make_pair(bid1, bid2)].second += 1;
                }
            }
        }

//    cout<<"updatedSCTrue size: "<<updatedSCTrue.size()<<" "<<updateSCSize<< endl;
//    cout<<"Initial overlayBatch size: "<<weightOverlay.size()<<endl;
        int sVertex;
        bool ifFound;
        for (auto it = updatedSCTrue.begin(); it != updatedSCTrue.end(); ++it) {
            bid1 = it->first.first, bid2 = it->first.second, wsum = it->second.first, cnt = it->second.second;
            ifFound = false;
            for (int i = 0; i < Tree[rank[bid1]].vert.size(); i++) {
                if (Tree[rank[bid1]].vert[i].first == bid2) {
                    ifFound = true;
                    if (Tree[rank[bid1]].vert[i].second.second < 1) {//need update
                        newdis = INF;
                        int countwt = 0;
                        for (auto it2 = Neighbor[bid1].begin(); it2 != Neighbor[bid1].end(); ++it2) {
                            if (it2->first == bid2) {
                                newdis = it2->second;//the weight value in the original graph
                                countwt = 1;
                                break;
                            }
                        }
                        int ssw, wtt, wid;
                        vector<pair<int, int>> Wnodes;
                        Wnodes.clear();
                        if (bid1 < bid2)
                            Wnodes = SCconNodesMT[bid1][bid2]; //cout<<"wid num "<<Wnodes.size()<<endl;
                        else
                            Wnodes = SCconNodesMT[bid2][bid1];
                        if (Wnodes.size() > 0) {
                            for (int i = 0; i < Wnodes.size(); i++) {
                                wid = Wnodes[i].first;
                                for (int j = 0; j < Tree[rank[wid]].vert.size(); j++) {
                                    if (Tree[rank[wid]].vert[j].first == bid1) {
                                        ssw = Tree[rank[wid]].vert[j].second.first;
                                    }
                                    if (Tree[rank[wid]].vert[j].first == bid2) {
                                        wtt = Tree[rank[wid]].vert[j].second.first;
                                    }
                                }

                                if (ssw + wtt < newdis) {
                                    newdis = ssw + wtt;
                                    countwt = 1;
                                } else if (ssw + wtt == newdis) {
                                    countwt += 1;
                                }
                            }
                        }

                        olddis = wsum;
                        if (newdis > olddis) {
//            sm->wait();
                            overlayBatch.emplace_back(make_pair(bid1, bid2), make_pair(olddis,
                                                                                       newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
                        } else if (newdis < olddis) {
//                            cout << "inc Something wrong happens. " << bid1 << " " << bid2 << " " << olddis << " " << newdis
//                                 << endl;//not for mix
//                            exit(1);
                        }
                    } else {
//                    cout<<"No update for this boundary shortcut! "<<bid1<<" "<<bid2<<" "<<Tree[rank[bid1]].vert[i].second.first<<" "<<Tree[rank[bid1]].vert[i].second.second<<" "<<cnt<<endl;
                    }
                    break;
                }
            }
            if (!ifFound) {
                cout << "Not found " << bid2 << " in " << bid1 << " 's vert. " << endl;
                exit(1);
            }
        }

    }

}

//Function for single-edge decrease update
void Graph::DecreaseSingle(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
            Neighbor[b][i].second=newW;
            break;
        }
    }
    for(int i=0;i<node_num;++i){
        vUpdated[i]=false;
    }
    map<int, vector<pair<pair<int,int>,pair<int,int>>>> partiBatch; partiBatch.clear();
    int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
    if(pid1!=pid2){//for cut edge
        DecreaseOverlay(a, b, newW,NeighborsOverlay,Tree,rank,heightMax);
    }else{//for intra edge
        vector<pair<pair<int,int>,pair<int,int>>> tempV;
        tempV.emplace_back(make_pair(a,b), make_pair(oldW,newW));
        partiBatch.insert({pid1, tempV});
        vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
        weightOverlay.clear();
        DecreaseParti(a,b,newW,NeighborsParti,Trees[pid1],ranks[pid1],heightMaxs[pid1]);

        //weightOverlay collect the changed edges on overlay graph
        vector<int> Bid=BoundVertex[pid1];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(int i=0;i<Bid.size();i++){
            bid1=Bid[i];
            for(int j=i+1;j<Bid.size();j++){
                bid2=Bid[j];
                if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                    olddis=NeighborsOverlay[bid1][bid2];
                }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                    continue;//exit(1);
                }

                newdis=QueryH2HPartition(bid1,bid2,pid1);
                if(newdis<olddis){
//                    cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
                    weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                }
            }
        }

        DecreaseOverlayBatch(weightOverlay,NeighborsOverlay,Tree,rank,heightMax,true);
        //cout<<"Overlay update number "<<weightOverlay.size()<<endl;
        //update the overlay graph index, after partition index update
        /*for(int l=0;l<weightOverlay.size();l++){
            Decrease(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay);
        }*/
    }
    // repair the partition index
    if(algoUpdate==PH2H_Post){
        Repair_PartiIndex(true, false, partiBatch);
    }

}

//Function for single-edge increase update
void Graph::IncreaseSingle(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
            Neighbor[b][i].second=newW;
            break;
        }
    }

    map<int, vector<pair<pair<int,int>,pair<int,int>>>> partiBatch; partiBatch.clear();

    for(int i=0;i<node_num;++i){
        vUpdated[i]=false;
    }

    int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
    //cout<<"increase edge in partition "<<pid<<endl;
    if(pid1!=pid2) {//for cut edge
        cout<<"Inter edge update"<<endl;
        IncreaseOverlay(a, b, oldW, newW,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid);
    }
    else {
        vector<pair<pair<int,int>,pair<int,int>>> tempV;
        tempV.emplace_back(make_pair(a,b), make_pair(oldW,newW));
        partiBatch.insert({pid1, tempV});
//        cout<<"zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"<<endl;
        vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
        weightOverlay.clear();

        IncreaseParti(a,b,oldW,newW,NeighborsParti,Trees[pid1],ranks[pid1],heightMaxs[pid1],SCconNodesMTP,VidtoTNidP);

        //cout<<"/////////////////////////////////////////"<<endl;

        //cout<<"boundary edge checkkkkkkkkkkkkkkkkkkkkkkkkkkk"<<endl;
        //boundary edges check
        vector<int> Bid=BoundVertex[pid1];
        int bid1,bid2,olddis,newdis;
        for(int i=0;i<Bid.size();i++){
            bid1=Bid[i];
            for(int j=i+1;j<Bid.size();j++){
                bid2=Bid[j];
                if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                    olddis=NeighborsOverlay[bid1][bid2];//only works for no-boundary
//                    olddis= QueryCore(bid1,bid2);
                }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                    continue;//exit(1);
                }
                newdis=QueryH2HPartition(bid1,bid2,pid1);
//                newdis=ShortcutDisCheck(bid1,bid2);
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                int overlaydis=QueryCore(bid1,bid2);
                if(newdis>olddis)//if '=', not problem; if '<', problem
                    weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                else if(newdis<olddis){
                    cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                    exit(1);
                }


            }
        }

        IncreaseOverlayBatch(weightOverlay,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);
//        CorrectnessCheckCore(100);
        //update the overlay graph index, after partition index update
        /*cout<<"Overlay update number "<<weightOverlay.size()<<endl;
        for(int l=0;l<weightOverlay.size();l++){
            Increase(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second.first,weightOverlay[l].second.second,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay,SCconNodesOverlayMT,VidtoTNidOverlay);
        }*/
        //cout<<"''''''''''''''''''''''''''''''''''''''''''"<<endl;
    }
    if(algoUpdate==PH2H_Post){
//        Trees=TreesNo;
        Repair_PartiIndex(true, true, partiBatch);
    }
}



