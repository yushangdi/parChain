#ifndef FINDNN_H
#define FINDNN_H

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include "shared.h"
#include "kdTree2.h"
#include "ndHash.h"
#include "bruteforce.h"
#include "node.h"

using namespace std;

namespace FINDNN {
    //find the farthest pair of points in P[start:start+n]

    /////// DualTree Traversal
    

    // computes shortest distance between two bounding boxes
    template<int dim, class nodeT>
    inline double node_distance(nodeT *Q, nodeT *R) {
        if (Q==R) return 0;
        for (int d = 0; d < dim; ++ d) {
        if (Q->pMin[d] > R->pMax[d] || R->pMin[d] > Q->pMax[d]) {
        // disjoint at dim d, and intersect on dim < d
        double rsqr = 0;
        for (int dd = d; dd < dim; ++ dd) {
        double tmp = max(Q->pMin[dd] - R->pMax[dd], R->pMin[dd] - Q->pMax[dd]);
        tmp = max(tmp, (double)0);
        rsqr += tmp * tmp;
        }
        return sqrt(rsqr);
        }
        }
        return 0; // intersect
    }

    template<int dim, class nodeT, class pointT>
    inline double node_distance(pointT *q, nodeT *R) {
        // if (Q==R) return 0;
        for (int d = 0; d < dim; ++ d) {
        if (q->p[d] > R->pMax[d] || R->pMin[d] > q->p[d]) {
        // disjoint at dim d, and intersect on dim < d
        double rsqr = 0;
        for (int dd = d; dd < dim; ++ dd) {
        double tmp = max(q->p[dd] - R->pMax[dd], R->pMin[dd] - q->p[dd]);
        tmp = max(tmp, (double)0);
        rsqr += tmp * tmp;
        }
        return sqrt(rsqr);
        }
        }
        return 0; // intersect
    }

    template<int dim, class nodeT>
    inline double node_far_distance(nodeT *Q, nodeT *R) {
        double result = 0;
        for (int d = 0; d < dim; ++ d) {
        double tmp = max(Q->pMax[d], R->pMax[d]) - min( Q->pMin[d], R->pMin[d]);
        result += tmp *tmp;
        }
        return sqrt(result);
    }

    // find the farthest points between two trees
    template<int dim, class nodeT>
    struct NNcomplete{
        LDS::EDGE *e;
        UnionFind::ParUF<intT> *uf;
        double ub = numeric_limits<double>::max();
        const LDS::edgeComparatorMax EC = LDS::edgeComparatorMax(); // use EC2?

        ~NNcomplete(){
            delete e;
        }

        NNcomplete(LDS::EDGE *ee, UnionFind::ParUF<intT> *t_uf):uf(t_uf),e(ee){}

        NNcomplete(LDS::EDGE *ee, UnionFind::ParUF<intT> *t_uf, double t_ub):uf(t_uf),e(ee), ub(t_ub){}


        NNcomplete(UnionFind::ParUF<intT> *t_uf):uf(t_uf){
            e = new LDS::EDGE(-1,-1,-1);
        }
        
        inline bool Score(double d, nodeT *Q, nodeT *R){
            return d < e->getW();
        }

        inline bool Score(nodeT *Q, nodeT *R, bool check = true){
            return check && Score(node_far_distance<dim, nodeT>(Q,R), Q, R);
        }

        inline void BaseCase(nodeT *Q, nodeT *R, intT i, intT j){
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            if(qrdist > ub){
                utils::writeMin(e, LDS::EDGE(-1,-1,numeric_limits<double>::max()), EC);
            }else{
                utils::writeMin(e, LDS::EDGE(Q->items[i]->idx(),R->items[j]->idx(),qrdist), EC);
            }
        }

        inline void BaseCaseSerial(nodeT *Q, nodeT *R, intT i, intT j){
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            if(qrdist > ub){
                e->update(-1,-1,numeric_limits<double>::max());
            }

            if(qrdist > e->getW()){
                e->update(Q->items[i]->idx(),R->items[j]->idx(),qrdist);
            }
        }

        inline double NodeDistForOrder(nodeT *Q, nodeT *R){
            return node_far_distance<dim, nodeT>(Q, R);
        }

        //if true, r first
        inline bool SpawnOrder(double l, double r){
            return l < r;
        }

        inline intT SpawnOrder(intT i){
            return 3-i;
        }

        inline void BasePost(nodeT *Q, nodeT *R){
            for(intT i=0; i<Q->size(); ++i){
                for(intT j=0;j<R->size() ; ++j){
                    BaseCase(Q,R,i,j);
                }
            }
        }
        inline void QLPost(nodeT *Q, nodeT *R){}
        inline void RLPost(nodeT *Q, nodeT *R){}
        inline void Post(nodeT *Q, nodeT *R){}
    };

    template<int dim, class nodeT>
    struct NNcomplete1{
        LDS::EDGE __attribute__ ((aligned (16))) e;
        double ub = numeric_limits<double>::max();
        const LDS::edgeComparatorMax EC = LDS::edgeComparatorMax();

        ~NNcomplete1(){
        }

        NNcomplete1(LDS::EDGE ee):e(ee){}

        NNcomplete1(LDS::EDGE ee, double t_ub):e(ee), ub(t_ub){}


        NNcomplete1(){
            e = LDS::EDGE(-1,-1,-1);
        }

        inline LDS::EDGE getResult(){
            return e;
        }

        inline double getResultW(){
            return e.getW();
        }

        inline bool isLeaf(nodeT *Q){
            return Q->isLeaf(); // || Q->size() < 400
        }
        
        inline bool Score(double d, nodeT *Q, nodeT *R){
            return d < e.getW();
        }

        inline bool Score(nodeT *Q, nodeT *R, bool check = true){
            return check && Score(node_far_distance<dim, nodeT>(Q,R), Q, R);
        }

        inline void BaseCase(nodeT *Q, nodeT *R, intT i, intT j){
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            if(qrdist > ub){
                utils::writeMin(&e, LDS::EDGE(-1,-1,numeric_limits<double>::max()), EC);
            }else{
                utils::writeMin(&e, LDS::EDGE(Q->items[i]->idx(),R->items[j]->idx(),qrdist), EC);
            }
        }


        inline double NodeDistForOrder(nodeT *Q, nodeT *R){
            return node_far_distance<dim, nodeT>(Q, R);
        }

        //if true, r first
        inline bool SpawnOrder(double l, double r){
            return l < r;
        }

        inline intT SpawnOrder(intT i){
            return 3-i;
        }

        inline void BasePost(nodeT *Q, nodeT *R){
            for(intT i=0; i<Q->size(); ++i){
                for(intT j=0;j<R->size() ; ++j){
                    BaseCase(Q,R,i,j);
                }
            }
        }
        inline void QLPost(nodeT *Q, nodeT *R){}
        inline void RLPost(nodeT *Q, nodeT *R){}
        inline void Post(nodeT *Q, nodeT *R){}
    };

    template<int dim, class nodeT, class objT>
    struct NNcomplete2{
        LDS::EDGE __attribute__ ((aligned (16))) e;
        double ub = numeric_limits<double>::max(); // used in naive computation
        const LDS::edgeComparatorMax EC = LDS::edgeComparatorMax(); 

        ~NNcomplete2(){
        }

        NNcomplete2(LDS::EDGE ee):e(ee){}

        NNcomplete2(LDS::EDGE ee, double t_ub):e(ee), ub(t_ub){}


        NNcomplete2(){
            e = LDS::EDGE(-1,-1,-1);
        }

        inline bool isLeaf(nodeT *Q){
            return Q->isLeaf(); // || Q->size() < 400
        }
        
        inline bool Score(double d, nodeT *Q, nodeT *R){
            return d < e.getW();
        }

        inline bool Score(nodeT *Q, nodeT *R, bool check = true){
            return check && Score(node_far_distance<dim, nodeT>(Q,R), Q, R);
        }

        inline double NodeDistForOrder(nodeT *Q, nodeT *R){
            return node_far_distance<dim, nodeT>(Q, R);
        }

        //if true, r first
        inline bool SpawnOrder(double l, double r){
            return l < r;
        }

        inline intT SpawnOrder(intT i){
            return 3-i;
        }

        inline void BasePost(nodeT *Q, nodeT *R){
            LDS::node<objT> *items1 = Q->items;
            LDS::EDGE locale = LDS::EDGE(-1,-1, -1);
            bool flag = true;

            while(true && flag){
                LDS::node<objT> *items2 = R->items;
                while(true && flag){
                    double qrdist = items1->elt->pointDist(items2->elt);
                    if(qrdist > ub){ // larger than upper bound, stop
                        utils::writeMin(&e, LDS::EDGE(-1,-1,numeric_limits<double>::max()), EC);
                        flag = false;
                        break;
                    }
                    if(qrdist > locale.getW()){
                        locale = LDS::EDGE(items1->elt->idx(),items2->elt->idx(),qrdist);
                    }
                    if(items2 == R->tail) break;
                    items2 = items2->next;  
                }
                if(items1 == Q->tail) break;
                items1 = items1->next;
            }
            utils::writeMin(&e, locale, EC);
        }

        inline void QLPost(nodeT *Q, nodeT *R){}
        inline void RLPost(nodeT *Q, nodeT *R){}
        inline void Post(nodeT *Q, nodeT *R){}
    };


    template<int dim, class nodeT>
    struct AllPtsNN{
        LDS::EDGE *edges;
        LDS::edgeComparator2 EC2;
#ifdef PERF_RANGE
        long *distance_computed;
        long *pointsInRange;
        unsigned long long *pointsInDist;

        void setCounter(long *t_distance_computed, long *t_pointsInRange, unsigned long long *t_pointsInDist){
            distance_computed = t_distance_computed;
            pointsInRange = t_pointsInRange;
            pointsInDist = t_pointsInDist;

        }
#endif
        AllPtsNN(LDS::EDGE *ee, double eps){
            EC2 = LDS::edgeComparator2(eps);
            edges = ee;
        }

        inline bool isLeaf(nodeT *Q){
            return Q->isLeaf();
        }

        inline bool Score(double d, nodeT *Q, nodeT *R){
            return d > (Q->nInfo).getUB();
        }

        inline bool Score(nodeT *Q, nodeT *R, bool check = true){
            return check && Score(node_distance<dim, nodeT>(Q,R), Q, R);
        }

        inline void BaseCase(nodeT *Q, nodeT *R, intT i, intT j){
#ifdef PERF_RANGE
            distance_computed[getWorkerId()*ELTPERCACHELINE]+=1;
            pointsInRange[getWorkerId()*ELTPERCACHELINE]+=1;
            pointsInDist[getWorkerId()*ELTPERCACHELINE]+=1;
#endif
            intT ii = Q->items[i]->idx();
            intT jj = R->items[j]->idx();
            if (ii == jj) return;
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            utils::writeMin(&edges[ii], LDS::EDGE(ii,jj,qrdist), EC2);
        }

        inline double NodeDistForOrder(nodeT *Q, nodeT *R){
            return node_distance<dim, nodeT>(Q, R);
        }

        //if true, r first
        inline bool SpawnOrder(double l, double r){
            return l > r;
        }

        inline intT SpawnOrder(intT i){
            return i;
        }

        void BasePost(nodeT *Q, nodeT *R){
            for(intT i=0; i<Q->size(); ++i){
                for(intT j=0;j<R->size() ; ++j){
                    BaseCase(Q,R,i,j);
                }
            }

            double temp = edges[Q->items[0]->idx()].getW();
            for(intT i=1; i<Q->size(); ++i){
               double eweight = edges[Q->items[i]->idx()].getW(); 
		    if(eweight > temp){
                    temp = eweight;
                }
            }
            (Q->nInfo).updateUB(temp);
        }

        inline void QLPost(nodeT *Q, nodeT *R){}
        inline void RLPost(nodeT *Q, nodeT *R){
            (Q->nInfo).updateUB(max((Q->left->nInfo).getUB(), (Q->right->nInfo).getUB()));
        }
        inline void Post(nodeT *Q, nodeT *R){
            RLPost(Q,R);
        }
    };

    template<int dim, class nodeT, class distF>
    struct AllPtsNNMetric{
        LDS::EDGE *edges;
        UnionFind::ParUF<intT> *uf;
        distF *distComputer;
        LDS::edgeComparator2 EC2;

        AllPtsNNMetric(LDS::EDGE *ee, UnionFind::ParUF<intT> *uuf, distF *ddistComputer, double eps):
            edges(ee), uf(uuf), distComputer(ddistComputer){
                EC2 = LDS::edgeComparator2(eps);
            }
        AllPtsNNMetric(double eps){
            EC2 = LDS::edgeComparator2(eps);
        }

        ~AllPtsNNMetric(){}


        inline intT pointIdxToEdgeIdx(nodeT *Q, intT i){
            return uf->find(Q->items[i]->idx());
        }

        inline bool isLeaf(nodeT *Q){
            return Q->isLeaf();
        }

        inline bool Score(double d, nodeT *Q, nodeT *R){
            return (R->nInfo.getCId() != -1 && R->nInfo.getCId() == Q->nInfo.getCId()) || d > (Q->nInfo).getUB();
        }

        inline bool Score(nodeT *Q, nodeT *R, bool check = true){
            return check && Score(node_distance<dim, nodeT>(Q,R), Q, R);
        }

        inline void BaseCase(nodeT *Q, nodeT *R, intT i, intT j){
            intT ii = Q->items[i]->idx();
            intT jj = R->items[j]->idx();
            intT cid1 = uf->find(ii);
            intT cid2 = uf->find(jj);
            if (cid1 == cid2) return;
            double qrdist = distComputer->getPointDist(Q->items[i], R->items[j]);//(Q->items[i])->pointDist(*(R->items[j]));
            utils::writeMin(&edges[cid1], LDS::EDGE(cid1,cid2,qrdist), EC2);
        }

        inline void BaseCase(nodeT *Q, nodeT *R, intT ii, intT jj, intT cid1){
            intT cid2 = uf->find(jj);
            if (cid1 == cid2) return;
            double qrdist = distComputer->getPointDist(ii, jj);//(Q->items[i])->pointDist(*(R->items[j]));
            utils::writeMin(&edges[cid1], LDS::EDGE(cid1,cid2,qrdist), EC2);
        }

        inline double NodeDistForOrder(nodeT *Q, nodeT *R){
            return node_distance<dim, nodeT>(Q, R);
        }

        //if true, r first
        inline bool SpawnOrder(double l, double r){
            return l > r;
        }

        inline intT SpawnOrder(intT i){
            return i;
        }

        void BasePost(nodeT *Q, nodeT *R){
            for(intT i=0; i<Q->size(); ++i){
                
                for(intT j=0;j<R->size() ; ++j){
                    BaseCase(Q,R,i,j);
                }
            }

            double temp = edges[pointIdxToEdgeIdx(Q, 0)].getW();
            for(intT i=1; i<Q->size(); ++i){
               double eweight = edges[pointIdxToEdgeIdx(Q, i)].getW(); 
		    if(eweight > temp){
                    temp = eweight;
                }
            }
            (Q->nInfo).updateUB(temp);
        }

        inline void QLPost(nodeT *Q, nodeT *R){}
        inline void RLPost(nodeT *Q, nodeT *R){
            (Q->nInfo).updateUB(max((Q->left->nInfo).getUB(), (Q->right->nInfo).getUB()));
        }
        inline void Post(nodeT *Q, nodeT *R){
            RLPost(Q,R);
        }
    };

    // find the NN of cluster cid
    //IMPORTANT: cid must be a valid cluster id, otherwise will find itself
    template<int dim, class nodeT>
    struct NNsingle{
        LDS::EDGE *e;
        UnionFind::ParUF<intT> *uf;
        intT cid; // cid of query cluster
        LDS::edgeComparator2 EC2;

        ~NNsingle(){
            delete e;
        }

        NNsingle(LDS::EDGE *t_e, UnionFind::ParUF<intT> *t_uf, intT t_cid, double eps):e(t_e),uf(t_uf), cid(t_cid){
            EC2 = LDS::edgeComparator2(eps);
        }

        NNsingle(UnionFind::ParUF<intT> *t_uf, intT t_cid, double eps):uf(t_uf), cid(t_cid){
            e = new LDS::EDGE(-1,-1,numeric_limits<double>::max());
            EC2 = LDS::edgeComparator2(eps);
        }

        inline bool isLeaf(nodeT *Q){
            return Q->isLeaf(); // || Q->size() < 400
        }
        
        inline bool Score(double d, nodeT *Q, nodeT *R){
            return (R->nInfo.getCId() == cid) || (d > e->getW());
        }

        inline bool Score(nodeT *Q, nodeT *R, bool check = true){
            return check && Score(node_distance<dim, nodeT>(Q,R), Q, R);
        }

        inline void BaseCase(nodeT *Q, nodeT *R, intT i, intT j){
            if((!Q->items[i] )|| (!R->items[j]) ) return;
            if(uf->find(R->items[j]->idx()) == cid) return;
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            utils::writeMin(e, LDS::EDGE(Q->items[i]->idx(),R->items[j]->idx(),qrdist), EC2);
        }

        inline void BaseCaseSerial(nodeT *Q, nodeT *R, intT i, intT j){
            if(uf->find(R->items[j]->idx()) == cid) return;
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            if(qrdist < e->getW()){
                e->update(Q->items[i]->idx(),R->items[j]->idx(),qrdist);
            }
            if(qrdist == e->getW() && uf->find(R->items[j]->idx()) < uf->find(e->second)){
                e->update(Q->items[i]->idx(),R->items[j]->idx(),qrdist);
            }
        }

        inline double NodeDistForOrder(nodeT *Q, nodeT *R){
            return node_distance<dim, nodeT>(Q, R);
        }

        //if true, r first
        inline bool SpawnOrder(double l, double r){
            return l > r;
        }

        inline intT SpawnOrder(intT i){
            return i;
        }

        inline void BasePost(nodeT *Q, nodeT *R){
            for(intT i=0; i<Q->size(); ++i){
                for(intT j=0;j<R->size(); ++j){
                    BaseCase(Q,R,i,j);
                }
            }
        }
        inline void QLPost(nodeT *Q, nodeT *R){}
        inline void RLPost(nodeT *Q, nodeT *R){}
        inline void Post(nodeT *Q, nodeT *R){}
    };


        // find the NN of cluster cid
    //IMPORTANT: cid must be a valid cluster id, otherwise will find itself
    template<int dim, class nodeT, class objT>
    struct NNaverage{

        NNaverage(){ }

        inline bool isLeaf(nodeT *Q){
            return Q->size() < 500; //todo: tune
        }

        inline tuple<double, bool> Score(nodeT *Q, nodeT *R, bool check = true){
            return make_tuple(0, false);
        }

        inline double BasePost(nodeT *Q, nodeT *R){
            LDS::node<objT> *items1 = Q->items;
            double d = 0;

            while(true){
                LDS::node<objT> *items2 = R->items;
                while(true){
                    d += items1->elt->pointDist(items2->elt);
                    if(items2 == R->tail) break;
                    items2 = items2->next;  
                }
                if(items1 == Q->tail) break;
                items1 = items1->next;
            }
            return d;

        }
        inline double QLPost(nodeT *Q, nodeT *R, double r1, double r2){return r1 + r2;}
        inline double RLPost(nodeT *Q, nodeT *R, double r1, double r2){return r1 + r2;}
        inline double Post(nodeT *Q, nodeT *R, double r1, double r2, double r3, double r4){
            return r1+r2+r3+r4;
        }
    };

    //Q and R are cluster trees
    template<class nodeT, class F>
    void dualtree_serial(nodeT *Q, nodeT *R, F *f, bool check = true){

        if(f->Score(Q,R,check)) return;

        if(f->isLeaf(Q) && f->isLeaf(R)){
            f->BasePost(Q,R);
        }else if (f->isLeaf(Q)){
           double dLeft = f->NodeDistForOrder(Q, R->left);
	       double dRight = f->NodeDistForOrder(Q, R->right);
           if(f->SpawnOrder(dLeft , dRight)){
              if(!f->Score(dRight, Q, R->right)) dualtree_serial<nodeT, F>(Q, R->right, f, false);
              if(!f->Score(dLeft, Q, R->left)) dualtree_serial<nodeT, F>(Q, R->left, f, false);
           }else{
              if(!f->Score(dLeft, Q, R->left)) dualtree_serial<nodeT, F>(Q, R->left, f, false);
              if(!f->Score(dRight, Q, R->right)) dualtree_serial<nodeT, F>(Q, R->right, f, false);
           }
           f->QLPost(Q,R);
        }else if(f->isLeaf(R)){
           double dLeft = f->NodeDistForOrder(Q->left, R);
	       double dRight = f->NodeDistForOrder(Q->right, R);
           if(f->SpawnOrder(dLeft , dRight)){
              if(!f->Score(dRight, Q->right, R)) dualtree_serial<nodeT, F>(Q->right, R, f, false);
              if(!f->Score(dLeft, Q->left, R)) dualtree_serial<nodeT, F>(Q->left, R, f, false);
           }else{
              if(!f->Score(dLeft, Q->left, R)) dualtree_serial<nodeT, F>(Q->left, R, f, false);
              if(!f->Score(dRight, Q->right, R)) dualtree_serial<nodeT, F>(Q->right, R, f, false);
           }
            f->RLPost(Q,R);
        }else{
            pair<double, pair<nodeT *, nodeT *>> callOrder[4];
            callOrder[0] = make_pair(f->NodeDistForOrder(Q->left, R->left), make_pair(Q->left, R->left));
            callOrder[1] = make_pair(f->NodeDistForOrder(Q->right, R->left), make_pair(Q->right, R->left));
            callOrder[2] = make_pair(f->NodeDistForOrder(Q->left, R->right), make_pair(Q->left, R->right));
            callOrder[3] = make_pair(f->NodeDistForOrder(Q->right, R->right), make_pair(Q->right, R->right));
            sort(callOrder, callOrder + 4);
            for (int cc = 0; cc < 4; ++ cc) {
                int c = f->SpawnOrder(cc);
                nodeT *QQ = callOrder[c].second.first;
                nodeT *RR = callOrder[c].second.second;
                if(!f->Score(callOrder[c].first, QQ, RR)) dualtree_serial<nodeT, F>(QQ, RR, f, false);
            }
            f->Post(Q,R);
        }
    }

        //Q and R are cluster trees
    template<class nodeT, class F, class E>
    E dualtree_serial2(nodeT *Q, nodeT *R, F *f, bool check = true){

        E result;
        bool stop;
        tie(result, stop) = f->Score(Q,R,check);
        if(stop) return result;

        if(f->isLeaf(Q) && f->isLeaf(R)){
            result = f->BasePost(Q,R);
        }else if (f->isLeaf(Q)){
            E result1 = dualtree_serial2<nodeT, F, E>(Q, R->left, f, false);
            E result2 = dualtree_serial2<nodeT, F, E>(Q, R->right, f, false);
            result = f->QLPost(Q,R, result1, result2);
        }else if(f->isLeaf(R)){
            E result1 = dualtree_serial2<nodeT, F, E>(Q->left, R, f, false);
            E result2 = dualtree_serial2<nodeT, F, E>(Q->right, R, f, false);
            result = f->RLPost(Q,R, result1, result2);
        }else{
            E result1 = dualtree_serial2<nodeT, F, E>(Q->left, R->left, f, false);
            E result2 = dualtree_serial2<nodeT, F, E>(Q->left, R->right, f, false);
            E result3 = dualtree_serial2<nodeT, F, E>(Q->right, R->right, f, false);
            E result4 = dualtree_serial2<nodeT, F, E>(Q->right, R->left, f, false);
	        
            result = f->Post(Q,R, result1, result2, result3, result4);
        }
        return result;
    }

    /////// SingleTree Traversal
    // dummy marker
    template<class Fr>
    struct DummyMarker{
        typedef typename Fr::nodeInfo nodeInfoT;
        typedef typename Fr::kdnodeT kdnodeT;
        typedef typename Fr::nodeT nodeT;
        typedef typename nodeInfoT::infoT infoT;
        infoT initVal = nodeInfoT().initInfoVal();
        DummyMarker(){}
        DummyMarker(UnionFind::ParUF<intT> *t_uf, nodeT *t_nodes, intT *t_rootIdx){}
        inline bool doMark(intT C, intT round){ return false;}
        inline bool isTopDown(infoT info){return false;}
        inline void TopDownNode(kdnodeT *Q, infoT info){}
        inline void BottomUpNode(kdnodeT *Q, infoT info){}
        inline void BaseCase(kdnodeT *Q, infoT info){}
        inline infoT SwitchMode(kdnodeT *Q, infoT info){return info;}
        inline bool Par(kdnodeT *Q){return false;}
        inline bool Stop(kdnodeT *Q, infoT info){return true;}
    };

    template<int dim, class Fr>
    struct MarkClusterId{
        typedef typename Fr::nodeInfo nodeInfoT;
        typedef typename Fr::nodeT nodeT;
        typedef typename Fr::kdnodeT kdnodeT;
        // bool isTopDown;
        typedef intT infoT;
        infoT initVal = nodeInfoT().initInfoVal();

        UnionFind::ParUF<intT> *uf;

        MarkClusterId(UnionFind::ParUF<intT> *t_uf, nodeT *t_nodes, intT *t_rootIdx):uf(t_uf){}
        MarkClusterId(){}


        inline bool doMark(intT C, intT round){ return round > 5 && C > 1;}

        inline bool isTopDown(intT id){ return id != -1;}

        inline void TopDownNode(kdnodeT *Q, intT id){
            if(!isTopDown(id)) return;
            Q->nInfo.setCId(id);
        }

        inline void BottomUpNode(kdnodeT *Q, intT id){
            if(isTopDown(id)) return;
            intT cidl = Q->left->nInfo.getCId();
            if(cidl != -1){
                intT cidr = Q->right->nInfo.getCId();
                if(cidl == cidr)Q->nInfo.setCId(cidl);
            }
        }

        inline void BaseCase(kdnodeT *Q, intT id){
            if(isTopDown(id)){
                Q->nInfo.setCId(id);
            }else{
                id = uf->find(Q->items[0]->idx());
                for(intT i=0; i<Q->size(); ++i){
                    if(uf->find(Q->items[i]->idx())!= id){
                        Q->nInfo.setCId(-1);
                        return;
                    }
                }
                Q->nInfo.setCId(id);
            }
        }

        inline intT SwitchMode(kdnodeT *Q, intT id){
            if(isTopDown(id)) return id;
            intT cid = Q->nInfo.getCId();
            if(cid == -1) return -1;
            return uf->find(cid); 
        }

        inline bool Par(kdnodeT *Q){return Q->size() > 2000;}

        inline bool Stop(kdnodeT *Q, intT id){
            intT cid = Q->nInfo.getCId();
            return cid != -1 && cid == id;
        }

    };

    template<int dim, class Fr>
    struct MarkClusterIdResetUpper{
        typedef typename Fr::nodeInfo nodeInfoT;
        typedef typename Fr::nodeT nodeT;
        typedef typename Fr::kdnodeT kdnodeT;
        // bool isTopDown;
        typedef intT infoT;
        infoT initVal = nodeInfoT().initInfoVal();

        UnionFind::ParUF<intT> *uf;

        MarkClusterIdResetUpper(UnionFind::ParUF<intT> *t_uf, nodeT *t_nodes, intT *t_rootIdx):uf(t_uf){}
        MarkClusterIdResetUpper(){}


        inline bool doMark(intT C, intT round){ return true;}

        inline bool isTopDown(intT id){ return id != -1;}

        inline void TopDownNode(kdnodeT *Q, intT id){
            if(!isTopDown(id)) return;
            Q->nInfo.setCId(id);
            Q->nInfo.resetUB();
        }

        inline void BottomUpNode(kdnodeT *Q, intT id){
            if(isTopDown(id)) return;
            intT cidl = Q->left->nInfo.getCId();
            if(cidl != -1){
                intT cidr = Q->right->nInfo.getCId();
                if(cidl == cidr)Q->nInfo.setCId(cidl);
            }
            Q->nInfo.resetUB();
        }

        inline void BaseCase(kdnodeT *Q, intT id){
            Q->nInfo.resetUB();
            if(isTopDown(id)){
                Q->nInfo.setCId(id);
            }else{
                id = uf->find(Q->items[0]->idx());
                for(intT i=0; i<Q->size(); ++i){
                    if(uf->find(Q->items[i]->idx())!= id){
                        Q->nInfo.setCId(-1);
                        return;
                    }
                }
                Q->nInfo.setCId(id);
            }
        }

        inline intT SwitchMode(kdnodeT *Q, intT id){
            if(isTopDown(id)) return id;
            intT cid = Q->nInfo.getCId();
            if(cid == -1) return -1;
            return uf->find(cid); 
        }

        inline bool Par(kdnodeT *Q){return Q->size() > 2000;}

        inline bool Stop(kdnodeT *Q, intT id){
            return false;
        }

    };


  // an optimization for AllNN
  // when no known upperbound is known, use diameter of node
  // only applicable when multiple clusters are in the node (cluster_id of node == -1)
    template<int dim, class Fr, class distT>
    struct MarkClusterIdResetUpperDiam{
        typedef typename Fr::nodeInfo nodeInfoT;
        typedef typename Fr::nodeT nodeT;
        typedef typename Fr::kdnodeT kdnodeT;
        // bool isTopDown;
        typedef intT infoT;
        infoT initVal = nodeInfoT().initInfoVal();

        UnionFind::ParUF<intT> *uf;
        distT *distComputer;

        MarkClusterIdResetUpperDiam(UnionFind::ParUF<intT> *t_uf, nodeT *t_nodes, 
                                    intT *t_rootIdx, distT *t_distComputer):uf(t_uf), distComputer(t_distComputer){}
        MarkClusterIdResetUpperDiam(){}

        inline bool doMark(intT C, intT round){ return true;}
        inline bool isTopDown(intT id){ return id != -1;}
        inline void TopDownNode(kdnodeT *Q, intT id){}

        inline void BottomUpNode(kdnodeT *Q, intT id){
            if(isTopDown(id)){
                Q->nInfo.setCId(id);
            }else{
                intT cidl = Q->left->nInfo.getCId();
                if(cidl != -1){
                    intT cidr = Q->right->nInfo.getCId();
                    if(cidl == cidr)Q->nInfo.setCId(cidl);
                }
            }
            Q->nInfo.resetUB(); 
            (Q->nInfo).updateUB(max((Q->left->nInfo).getUB(), (Q->right->nInfo).getUB()));
            distComputer->updateUB(Q);
        }

        inline void BaseCase(kdnodeT *Q, intT id){ 
            if(isTopDown(id)){
                Q->nInfo.setCId(id);
            }else{
                id = uf->find(Q->items[0]->idx());
                for(intT i=0; i<Q->size(); ++i){
                    if(uf->find(Q->items[i]->idx())!= id){
                        id = -1;//Q->nInfo.setCId(-1);
                        break;
                    }
                }
                Q->nInfo.setCId(id);
            }
            Q->nInfo.resetUB(); 
            distComputer->updateUB(Q);
        }

        inline intT SwitchMode(kdnodeT *Q, intT id){
            if(isTopDown(id)) return uf->find(id);
            intT cid = Q->nInfo.getCId();
            if(cid == -1) return -1;
            return uf->find(cid); 
        }

        inline bool Par(kdnodeT *Q){return Q->size() > 2000;}

        inline bool Stop(kdnodeT *Q, intT id){
            return false;
        }

    };


  
    struct DummyMarkerCenters{
        DummyMarkerCenters(intT *t_sizes){
        }
        DummyMarkerCenters(){
        }
    };

    // mark cluster id and others on the tree of Centers
    template<int dim, class pointT, class nodeInfo>
    struct MarkKdTreeCenters{
        typedef typename nodeInfo::infoT infoT;
        typedef kdNode<dim, pointT, nodeInfo> kdnodeT;

        intT *sizes;
        infoT initVal = nodeInfo().initInfoVal();

        MarkKdTreeCenters(intT *t_sizes):sizes(t_sizes){
        }
        MarkKdTreeCenters(){
        }

        inline bool doMark(intT C, intT round){ return true;}

        inline bool isTopDown(infoT info){return false;}

        inline void TopDownNode(kdnodeT *Q, infoT info){
            if(!isTopDown(info)) return;
            Q->nInfo.setInfo(info);
        }

        inline void BottomUpNode(kdnodeT *Q, infoT info){
            if(isTopDown(info)) return;
            Q->nInfo.setMinN(min(Q->left->nInfo.getMinN(), Q->right->nInfo.getMinN()));
        }

        inline void BaseCase(kdnodeT *Q, infoT info){
            if(isTopDown(info)){
                Q->nInfo.setInfo(info);
            }else{
                intT id = Q->items[0]->idx();
                intT min_n = sizes[id]; 
                for(intT i=1; i<Q->size(); ++i){
                    intT id_temp = Q->items[i]->idx();
                    min_n = min(min_n, sizes[id_temp]);  
                }
                id = -1;
                Q->nInfo.setInfo(infoT(id, min_n));
            }
        }

        inline infoT SwitchMode(kdnodeT *Q, infoT info){ //must be -1
            return info;
        }

        inline bool Par(kdnodeT *Q){return Q->size() > 1000;}

        inline bool Stop(kdnodeT *Q, infoT info){
            return false;
        }

    };

    // mark cluster core distances
    template<int dim, class kdnodeT>
    struct MarkKdTreeCoreDist{
        typedef int infoT;

        double *coredists;
        infoT initVal = 0;//dummy

        MarkKdTreeCoreDist(double *t_coredists):coredists(t_coredists){
        }
        MarkKdTreeCoreDist(){
        }

        inline bool isTopDown(infoT info){return false;}//{ return get<0>(info) != -1;}

        inline void TopDownNode(kdnodeT *Q, infoT info){
        }

        inline void BottomUpNode(kdnodeT *Q, infoT info){
            Q->nInfo.max_core_dist = max((Q->left->nInfo).max_core_dist, (Q->right->nInfo).max_core_dist);
        }

        inline void BaseCase(kdnodeT *Q, infoT info){
            double temp = coredists[Q->items[0]->idx()];
            for(intT i=1; i<Q->size(); ++i){
            double q_cd = coredists[Q->items[i]->idx()]; 
                if(q_cd > temp) temp = q_cd;
            }
            Q->nInfo.max_core_dist = temp;
        }

        inline infoT SwitchMode(kdnodeT *Q, infoT info){ //must be -1
            return info;
        }

        inline bool Par(kdnodeT *Q){return Q->size() > 1000;}

        inline bool Stop(kdnodeT *Q, infoT info){
            return false;
        }
    };



    template<class nodeT, class F, class E>
    void singletree(nodeT *Q, F *f, E info){
        info = f->SwitchMode(Q, info);
        if(f->Stop(Q, info)) return;
        if(Q->isLeaf()){
            f->BaseCase(Q, info);
        }else{
            f->TopDownNode(Q, info);
            if(f->Par(Q)){
                cilk_spawn singletree<nodeT, F, E>(Q->left, f, info);
                singletree<nodeT, F, E>(Q->right, f, info);
                cilk_sync;
            }else{
                singletree<nodeT, F, E>(Q->left, f, info);
                singletree<nodeT, F, E>(Q->right, f, info);
            }
            f->BottomUpNode(Q, info);
        }
    }


    template<int dim, class nodeT>
    void printTreeCId(nodeT *Q){
        typedef UTIL::PrintClusterId<dim, nodeT> M;
        M *printer = new M();
        singletree<nodeT, M, intT>(Q, printer, -1);
    }
}

#endif