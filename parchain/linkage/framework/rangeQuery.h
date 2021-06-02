#ifndef CLINK_RANGE_QUERY_H
#define CLINK_RANGE_QUERY_H

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include "shared.h"
#include "kdTree2.h"
#include "ndHash.h"
#include "bruteforce.h"
#include "neighbor.h"
#include "neighbor_parallel.h"
#include "serialHash.h"
#include "box.h"

using namespace std;

namespace FINDNN {

    /////// Range Query
    template<int dim, class pointTT, class nodeInfoT, class nodeTT>
    struct DummyRangeF{
        typedef pointTT pointT;
        typedef nodeInfoT nodeInfo;
        typedef nodeTT nodeT;
        typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
    };
    

    // used for kdtree
    template<intT dim, class objT, class nodeInfoT, class distT, class Box>
    struct RangeQueryCountF1{
        typedef typename distT::pointT pointT;
        typedef nodeInfoT nodeInfo;
        typedef typename distT::nodeT nodeT;
        typedef typename distT::clusterCacheT clusterCacheT;
        typedef kdTree<dim, pointT, nodeInfo> kdtreeT;
        typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
    
        UnionFind::ParUF<intT> *uf;
        intT cid;
        pair<intT, double> e;
        clusterCacheT *tb; //
        intT pid;
        Box box;
        distT *distComputer;
        LDS::edgeComparator2 EC2;
        double eps;
        const bool local = true;

        RangeQueryCountF1(UnionFind::ParUF<intT> *t_uf, intT t_cid, 
            nodeT *t_nodes, intT *t_rootIdx, LDS::distCacheT **t_tbs, LDS::EDGE *t_edges,
            distT *t_distComputer, bool t_no_cache, intT C, double _eps):
            uf(t_uf), cid(t_cid), //edges(t_edges), 
            distComputer(t_distComputer), eps(_eps){
            EC2 = LDS::edgeComparator2(eps);
            e = make_pair(t_edges[cid].second, t_edges[cid].getW());

            pid = getWorkerId();
            tb = distComputer->initClusterTb(pid, C);//clusterTbs[idx];
            box = Box();
        }

        ~RangeQueryCountF1(){
        }

        inline intT getFinalNN(){return e.first;}
        inline double getFinalDist(){return e.second;}

        inline void updateDist(intT Rid, bool reach_thresh){
            if(cid != Rid && Rid != e.first){
                double dist = distComputer->getDistNaive(cid,Rid, -1, e.second, false); //, false

                if(e.second - dist > eps){ e = make_pair(Rid, dist);}  
                else if(abs(e.second - dist) <= eps && Rid < e.first){e = make_pair(Rid, dist); }
                // tb->deleteVal(Rid); //does not support delete and insert at the same time, need to remove if parallel
            }
        }

        inline tuple<intT, bool> incrementTable(intT Rid, intT a = 1){
            return distComputer->incrementTable(tb, Rid,  cid, a);
            // bool inserted; bool reach_thresh;
            // tie(inserted, reach_thresh) = tb->insert_thresh(make_pair(Rid,a));
            // if (inserted) return make_tuple(a, false);
            // if (reach_thresh) return make_tuple(0, true);
            // return make_tuple(tb->find_thresh(Rid).second, false);
        }

        inline bool isComplete(){return false;}
        inline bool isComplete2(kdnodeT *Q){
            intT  Rid = Q->nInfo.getCId();
            if(cid == Rid ) return true;
            if( Rid != -1){
                intT ct; bool reach_thresh;
                tie(ct, reach_thresh) = incrementTable(Rid, Q->size());
                if (reach_thresh || ct ==  distComputer->kdtrees[Rid]->getN()) updateDist(Rid, reach_thresh);
                return true;
            }else{
                return false;
            }
        }

        inline bool checkComplete(objT *p){
            // if(p->pointDist(qnode->center) > r + EC2.eps) return false;
            intT  Rid = uf->find(p->idx());
            if(cid == Rid ) return false;
            intT ct; bool reach_thresh;
            tie(ct, reach_thresh) = incrementTable(Rid);
            if (reach_thresh || ct ==  distComputer->kdtrees[Rid]->getN()) updateDist(Rid, reach_thresh);
            return false;
        }

        inline bool Par(kdnodeT *Q){
            return false;  // have to be false if using hashtable for clsuterhash
        }

        inline double getBall(nodeT* query, double r){
            return box.getBall(query, r);
        }

    };

    // need t_m active hash table size to store candidates
    // invariant: e contain the current nearest neighbor in tbs to cid
    // insert into (smallid, large id) table, only succeeded one compute and update
    template<intT dim, class objT, class nodeInfoT, class distT, class Box>
    struct RangeQueryCenterF{
        typedef typename distT::pointT pointT;
        typedef nodeInfoT nodeInfo;
        typedef typename distT::nodeT nodeT;
        typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
        typedef kdTree<dim, pointT, nodeInfo> kdtreeT;

        UnionFind::ParUF<intT> *uf;
        intT cid;
        LDS::EDGE *edges;
        nodeT *nodes;
        nodeT *qnode;
        intT *rootIdx;
        LDS::distCacheT **tbs; //
        LDS::edgeComparator2 EC2;
        distT *distComputer;
        Box box;
        double r;
        bool no_cache;
        const bool local = false; // writemin after


        RangeQueryCenterF(UnionFind::ParUF<intT> *t_uf, intT t_cid, 
            nodeT *t_nodes, intT *t_rootIdx, LDS::distCacheT **t_tbs, LDS::EDGE *t_edges,
            distT *t_distComputer, bool t_no_cache, intT C, double eps):
            uf(t_uf), cid(t_cid), //clusteredPts(t_clusteredPts),
            distComputer(t_distComputer),
            no_cache(t_no_cache){
                EC2 = LDS::edgeComparator2(eps);
            // keep nn candidate when merging
            tbs = t_tbs;
            edges = t_edges;
            nodes = t_nodes;
            rootIdx = t_rootIdx;
            // f = new F();
            qnode = getNode(cid);
            box = Box();
        }

        ~RangeQueryCenterF(){
        }
        inline intT getFinalNN(){return edges[cid].second;}
        inline double getFinalDist(){return edges[cid].getW();}

        inline intT idx(nodeT* node){ return node->idx;}
        inline nodeT *getNode(intT cid){return nodes+rootIdx[cid];}
        inline intT idx(intT cid){return idx(getNode(cid));}

        inline double my_node_distance_sq(kdnodeT *Q) {
            pointT qcenter = qnode->center;
            for (int d = 0; d < dim; ++ d) {
                if (Q->pMin[d] > qcenter[d] || qcenter[d] > Q->pMax[d]) {
                // disjoint at dim d, and intersect on dim < d
                double rsqr = 0;
                for (int dd = d; dd < dim; ++ dd) {
                    double tmp = max(Q->pMin[dd] - qcenter[dd], qcenter[dd] - Q->pMax[dd]);
                    tmp = max(tmp, (double)0);
                    rsqr += tmp * tmp;
                }
                return rsqr;
                }
            }
            return 0; // intersect
        }


        // return 0 if not found
        // return distance if found
        inline double find(intT qid, intT rid){
            CHECK_NO_CACHE(-203)
            intT qIdx = idx(qid);
            intT rIdx = idx(rid);
            
            typename LDS::distCacheT::eType result;
            bool reach_thresh;
            tie(result, reach_thresh) = tbs[qIdx]->find_thresh(rid);
            if(!reach_thresh && result.idx == rIdx){
			return result.dist;
            }
            
            tie(result, reach_thresh) = tbs[rIdx]->find_thresh(qid);
            if(!reach_thresh && result.idx == qIdx){
			return result.dist;
            }
            return UNFOUND_TOKEN;
        }

        inline void insert(intT qid, intT rid, double d){
            if(d == LARGER_THAN_UB){
                return;
            }
            CHECK_NO_CACHE(-222)
            intT qIdx = idx(qid);
            intT rIdx = idx(rid);

            tbs[qIdx]->insert2(LDS::hashClusterAveET(rid, rIdx, d));
            tbs[rIdx]->insert2(LDS::hashClusterAveET(qid, qIdx, d));
        }

        // return true when insert if sucussful or tables full
        // return true means need to compute distance
        inline bool insert_check(intT qid, intT rid){
            CHECK_NO_CACHE(-233)
            if(qid > rid){
                swap(qid, rid);
            }
            intT qIdx = idx(qid);
            intT rIdx = idx(rid);
            bool inserted; bool reach_thresh;
            tie(inserted, reach_thresh) = tbs[qIdx]->insert_thresh(LDS::hashClusterAveET(rid, rIdx, CHECK_TOKEN));
            if(!reach_thresh) return inserted;

            tie(inserted, reach_thresh) = tbs[rIdx]->insert_thresh(LDS::hashClusterAveET(qid, qIdx, CHECK_TOKEN));
            if(!reach_thresh) return inserted;
            return true;
        }

        inline void updateDist(intT Rid){ // need another table!
            // only first inserted is true, CHECK_TOKEN as a placeholder
            // if already in table, entry will not be replaced because CHECK_TOKEN is inserted
            // success = true when pair not in tbs and this is the first insert
                // if((cid == 609948 && Rid == 563148) || (cid ==  563148 && Rid == 609948)){
                //     cout << cid << " !!!!!" <<endl;
                // }
            if(!no_cache){
            bool success = insert_check(cid, Rid);
            if(!success){  // only compute distance once
                double dist = find(cid, Rid);
                // success = false only when insertions fail and reach_thresh is false
                if(dist == UNFOUND_TOKEN){  cout << "should not find unfound_token" << endl;exit(1);}
                if(dist != CHECK_TOKEN){              
                    utils::writeMin(&edges[cid], LDS::EDGE(cid, Rid, dist), EC2);
                    utils::writeMin(&edges[Rid], LDS::EDGE(Rid, cid, dist), EC2);
                }
                return;  
            }
            }
            
            double dist;
            if(distComputer->id_only) {dist = distComputer->getDistNaive(cid, Rid);}
            else {dist = distComputer->getDistNaive(qnode, getNode(Rid));}
            if(!no_cache) insert(cid, Rid, dist); //TODO: if both full, do not insert
            //in case Rid searches for cid
            utils::writeMin(&edges[cid], LDS::EDGE(cid, Rid, dist), EC2);
            utils::writeMin(&edges[Rid], LDS::EDGE(Rid, cid, dist), EC2);

        }

        inline bool isComplete(){return false;}
        inline bool isComplete2(kdnodeT *Q){
            if(distComputer->method == LDS::WARD){
                double dsq = my_node_distance_sq(Q);
                double min_n = (double)Q->nInfo.getMinN();
                double qn = (double)qnode->size();
                if(dsq > (qn + min_n)/min_n / qn  / 2 * r * r ) return true;
            }
            return false;
        }

        inline bool checkComplete(objT *p){
            //already checked by iteminball
            // if(p->pointDist(qnode->center) > r) return false; //eps already added in get box
            intT  Rid = p->idx(); //should all be centers uf->find(p->idx());
            if(cid != Rid && Rid != edges[cid].second) updateDist(Rid);
            return false;
        }

        inline bool Par(kdnodeT *Q){
            return Q->size() > 2000; //TODO: try turn off
            // return false; 
        }

        // find points in `root` within box of side `r` centered at `query` node
        inline double getBall(nodeT* query, double _r){
            r = _r;
            return box.getBall(query, r);
        }

    };

    template<class pointT, class func, class nodeT>
    void NNcandidate(nodeT *Q, pointT center, double r, func* f) {
    int relation = Q->boxBallCompare(center, r, Q->pMin, Q->pMax);
    if(relation == Q->boxExclude) return;
    if (f->isComplete2(Q)) return;

    if (Q->isLeaf()) {
        for(intT i=0; i<Q->size(); ++i) {
            if (Q->itemInBall(center, r, Q->items[i])) {
                f->checkComplete(Q->items[i]);
            }
        }
    } else {
        if(f->Par(Q)){
            cilk_spawn NNcandidate(Q->left, center, r, f);
            NNcandidate(Q->right, center, r, f);
            cilk_sync;
        }else{
            NNcandidate(Q->left, center, r, f);
            NNcandidate(Q->right, center, r, f);
        }
    }
  }

    // specific isLeaf
    template<int dim, class nodeT, class F>
    void dualtree2_serial(nodeT *Q, nodeT *R, F *f, bool check = true){
        if(f->Score(Q,R,check)) return;

        if(f->isLeafQ(Q) && f->isLeafR(R)){
            // for(intT i=0; i<Q->size(); ++i){
            //     for(intT j=0;j<R->n; ++j){
            //         f->BaseCase(Q,R,i,j);
            //     }
            // }
            f->BasePost(Q,R);
        }else if (f->isLeafQ(Q)){
            dualtree2_serial<dim, nodeT, F>(Q, R->right, f, false);
            dualtree2_serial<dim, nodeT, F>(Q, R->left, f, false);
        }else if(f->isLeafR(R)){
            dualtree2_serial<dim, nodeT, F>(Q->right, R, f, false);
            dualtree2_serial<dim, nodeT, F>(Q->left, R, f, false);
        }else{
            pair<double, pair<nodeT *, nodeT *>> callOrder[4];
            callOrder[0] = make_pair(f->NodeDistForOrder(Q->left, R->left), make_pair(Q->left, R->left));
            callOrder[1] = make_pair(f->NodeDistForOrder(Q->right, R->left), make_pair(Q->right, R->left));
            callOrder[2] = make_pair(f->NodeDistForOrder(Q->left, R->right), make_pair(Q->left, R->right));
            callOrder[3] = make_pair(f->NodeDistForOrder(Q->right, R->right), make_pair(Q->right, R->right));
            // sort(callOrder, callOrder + 4);
            for (int cc = 0; cc < 4; ++ cc) {
                int c = f->SpawnOrder(cc);
                nodeT *QQ = callOrder[c].second.first;
                nodeT *RR = callOrder[c].second.second;
                dualtree2_serial<dim, nodeT, F>(QQ, RR, f, false);
                // if(!f->Score(callOrder[c].first, QQ, RR)) dualtree2_serial<dim, nodeT, F>(QQ, RR, f, false);
            }
            f->Post(Q,R);
        }
    }

}
 
#endif
