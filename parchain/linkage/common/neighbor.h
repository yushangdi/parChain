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

using namespace std;

namespace FINDNN {
    //find the farthest pair of points in P[start:start+n]


    struct CLinkNodeInfo{
        intT cId;
        // double lb;
        double ub;
        intT round = 0;
        intT idx = -1; // node idx

        CLinkNodeInfo(){
            cId = -1;
            // lb = numeric_limits<double>::lowest();
            ub = numeric_limits<double>::max(); // used for allptsnn
            // round = 0;
        }

        inline double getUB(){return ub;}
        inline void updateUB(double tmp){
            utils::writeMin(&ub, tmp);
        }
        inline intT getCId(){return cId;}
        inline void setCId(intT id){cId = id;}
        inline intT getRound(){return round;}
        inline void setRound(intT r){round = r;}
        inline intT getIdx(){return idx;}
        inline void setIdx(intT r){idx = r;}
    };

    template<int dim, class objT>
    struct AveLinkNodeInfo2{
        typedef point<dim> pointT;
        intT cId;
        intT round = 0;
        intT idx; // node idx
        intT n = 1;
        LDS::node<objT> *items;
        AveLinkNodeInfo2<dim, objT> *left = nullptr;
        AveLinkNodeInfo2<dim, objT> *right = nullptr;
        LDS::node<objT> *tail;
        pointT pMin, pMax;

        AveLinkNodeInfo2(intT  t_cid, intT t_round, intT t_idx, AveLinkNodeInfo2<dim, objT> *t_left, AveLinkNodeInfo2<dim, objT> *t_right):
            cId(t_cid),
            round(t_round),
            idx(t_idx),
            items(t_left->items),
            left(t_left),
            right(t_right){
                n = t_left->n + t_right->n;
                tail = t_right->tail;
                left->tail->next = t_right->items;
                pMin = pointT(left->pMin.x);
                pMax = pointT(left->pMax.x);
                pMin.minCoords(right->pMin);
                pMax.maxCoords(right->pMax);
            }

        AveLinkNodeInfo2(intT  t_cid, LDS::node<objT> *t_items):
            cId(t_cid),
            idx(t_cid),
            items(t_items),
            tail(t_items){
                pMin = t_items->elt->p;
                pMax = t_items->elt->p;
            }

        inline bool isLeaf(){return n == 1;}
    };

    template<int dim, class objT>
    struct AveLinkNodeInfo1{
        intT cId;
        intT round = 0;
        intT idx; // node idx
        intT n = 1;
        AveLinkNodeInfo1<dim, objT> *left = nullptr;
        AveLinkNodeInfo1<dim, objT> *right = nullptr;
        point<dim> pMin, pMax;
        intT offset = -1; 

        AveLinkNodeInfo1(intT  t_cid, intT t_round, intT t_idx, AveLinkNodeInfo1<dim, objT> *t_left, AveLinkNodeInfo1<dim, objT> *t_right):
            cId(t_cid),
            round(t_round),
            idx(t_idx),
            left(t_left),
            right(t_right){
                n = t_left->n + t_right->n;
                pMin = left->pMin;//pointT(left->pMin.x);
                pMax = left->pMax;//pointT(left->pMax.x);
                pMin.minCoords(right->pMin);
                pMax.maxCoords(right->pMax);
            }
        
        AveLinkNodeInfo1(intT  t_cid, iPoint<dim> p):
            cId(t_cid),
            idx(t_cid),
            offset(t_cid){
                pMin = p.p;
                pMax = p.p;
            }

        inline bool isLeaf(){return n == 1;}
        inline intT getOffset(){return offset;}
        inline void setOffset(intT i){offset = i;}

    };

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
                // utils::writeMin(e, LDS::EDGE(-1,-1,numeric_limits<double>::max()), EC);
                e->update(-1,-1,numeric_limits<double>::max());
            }
            // else{
            //     utils::writeMin(e, LDS::EDGE(Q->items[i]->idx(),R->items[j]->idx(),qrdist), EC);
            // }
            if(qrdist > e->getW()){
                e->update(Q->items[i]->idx(),R->items[j]->idx(),qrdist);
            }
            // if(qrdist == e->getW() && uf->find(R->items[j]->idx()) < uf->find(e->second)){
            //     e->update(Q->items[i]->idx(),R->items[j]->idx(),qrdist);
            // }
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

        inline void BasePost(nodeT *Q, nodeT *R){}
        inline void QLPost(nodeT *Q, nodeT *R){}
        inline void RLPost(nodeT *Q, nodeT *R){}
        inline void Post(nodeT *Q, nodeT *R){}
    };

    template<int dim, class nodeT>
    struct AllPtsNN{
        LDS::EDGE *edges;
        const LDS::edgeComparator2 EC2 = LDS::edgeComparator2();

        AllPtsNN(LDS::EDGE *ee){
            edges = ee;
        }

        inline bool Score(double d, nodeT *Q, nodeT *R){
            return d > (Q->nInfo).getUB();
        }

        inline bool Score(nodeT *Q, nodeT *R, bool check = true){
            return check && Score(node_distance<dim, nodeT>(Q,R), Q, R);
        }

        inline void BaseCase(nodeT *Q, nodeT *R, intT i, intT j){
            intT ii = Q->items[i]->idx();
            intT jj = R->items[j]->idx();
            if (ii == jj) return;
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            utils::writeMin(&edges[ii], LDS::EDGE(ii,jj,qrdist), EC2);
        }

        inline void BaseCaseSerial(nodeT *Q, nodeT *R, intT i, intT j){
            intT ii = Q->items[i]->idx();
            intT jj = R->items[j]->idx();
            if (ii == jj) return;
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            // utils::writeMin(&edges[ii], LDS::EDGE(ii,jj,qrdist), EC2);
            if(qrdist < edges[ii].getW()){
                edges[ii].update(ii,jj,qrdist);
            }
            if(qrdist == edges[ii].getW() && jj < edges[ii].second){
                edges[ii].update(ii, jj,qrdist);
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

        void BasePost(nodeT *Q, nodeT *R){
            double temp = edges[Q->items[0]->idx()].getW();
            for(intT i=1; i<Q->n; ++i){
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

    // find the NN of cluster cid
    //IMPORTANT: cid must be a valid cluster id, otherwise will find itself
    template<int dim, class nodeT>
    struct NNsingle{
        LDS::EDGE *e;
        UnionFind::ParUF<intT> *uf;
        intT cid; // cid of query cluster
        const LDS::edgeComparator2 EC2 = LDS::edgeComparator2();
        
        ~NNsingle(){
            delete e;
        }

        NNsingle(LDS::EDGE *t_e, UnionFind::ParUF<intT> *t_uf, intT t_cid):e(t_e),uf(t_uf), cid(t_cid){}

        NNsingle(UnionFind::ParUF<intT> *t_uf, intT t_cid):uf(t_uf), cid(t_cid){
            e = new LDS::EDGE(-1,-1,numeric_limits<double>::max());
        }
        
        inline bool Score(double d, nodeT *Q, nodeT *R){
            return (R->nInfo.getCId() == cid) || (d > e->getW());
        }

        inline bool Score(nodeT *Q, nodeT *R, bool check = true){
            return check && Score(node_distance<dim, nodeT>(Q,R), Q, R);
        }

        inline void BaseCase(nodeT *Q, nodeT *R, intT i, intT j){
            if(uf->find(R->items[j]->idx()) == cid) return;
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            utils::writeMin(e, LDS::EDGE(Q->items[i]->idx(),R->items[j]->idx(),qrdist), EC2);
        }

        inline void BaseCaseSerial(nodeT *Q, nodeT *R, intT i, intT j){
            if(uf->find(R->items[j]->idx()) == cid) return;
            double qrdist = (Q->items[i])->pointDist(*(R->items[j]));
            // utils::writeMin(e, LDS::EDGE(Q->items[i]->idx(),R->items[j]->idx(),qrdist), EC2);
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

        inline void BasePost(nodeT *Q, nodeT *R){}
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
            return Q->n < 500;
            // return Q->n == 1;
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
    template<int dim, class nodeT, class F>
    void dualtree_serial(nodeT *Q, nodeT *R, F *f, bool check = true){

        if(f->Score(Q,R,check)) return;

        if(Q->isLeaf() && R->isLeaf()){
            for(intT i=0; i<Q->n; ++i){
                for(intT j=0;j<R->n; ++j){
                    f->BaseCase(Q,R,i,j);
                }
            }
            f->BasePost(Q,R);
        }else if (Q->isLeaf()){
           double dLeft = f->NodeDistForOrder(Q, R->left);
	       double dRight = f->NodeDistForOrder(Q, R->right);
           if(f->SpawnOrder(dLeft , dRight)){
              if(!f->Score(dRight, Q, R->right)) dualtree_serial<dim, nodeT, F>(Q, R->right, f, false);
              if(!f->Score(dLeft, Q, R->left)) dualtree_serial<dim, nodeT, F>(Q, R->left, f, false);
           }else{
              if(!f->Score(dLeft, Q, R->left)) dualtree_serial<dim, nodeT, F>(Q, R->left, f, false);
              if(!f->Score(dRight, Q, R->right)) dualtree_serial<dim, nodeT, F>(Q, R->right, f, false);
           }
           f->QLPost(Q,R);
        }else if(R->isLeaf()){
           double dLeft = f->NodeDistForOrder(Q->left, R);
	       double dRight = f->NodeDistForOrder(Q->right, R);
           if(f->SpawnOrder(dLeft , dRight)){
              if(!f->Score(dRight, Q->right, R)) dualtree_serial<dim, nodeT, F>(Q->right, R, f, false);
              if(!f->Score(dLeft, Q->left, R)) dualtree_serial<dim, nodeT, F>(Q->left, R, f, false);
           }else{
              if(!f->Score(dLeft, Q->left, R)) dualtree_serial<dim, nodeT, F>(Q->left, R, f, false);
              if(!f->Score(dRight, Q->right, R)) dualtree_serial<dim, nodeT, F>(Q->right, R, f, false);
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
                if(!f->Score(callOrder[c].first, QQ, RR)) dualtree_serial<dim, nodeT, F>(QQ, RR, f, false);
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

    template<int dim, class nodeT>
    struct MarkClusterId{
        // bool isTopDown;
        UnionFind::ParUF<intT> *uf;

        MarkClusterId(UnionFind::ParUF<intT> *t_uf):uf(t_uf){}

        inline bool isTopDown(intT id){ return id != -1;}

        inline void TopDownNode(nodeT *Q, intT id){
            if(!isTopDown(id)) return;
            Q->nInfo.setCId(id);
        }

        inline void BottomUpNode(nodeT *Q, intT id){
            if(isTopDown(id)) return;
            intT cidl = Q->left->nInfo.getCId();
            if(cidl != -1){
                intT cidr = Q->right->nInfo.getCId();
                if(cidl == cidr)Q->nInfo.setCId(cidl);
            }
        }

        inline void BaseCase(nodeT *Q, intT id){
            if(isTopDown(id)){
                Q->nInfo.setCId(id);
            }else{
                id = uf->find(Q->items[0]->idx());
                for(intT i=0; i<Q->n; ++i){
                    if(uf->find(Q->items[i]->idx())!= id){
                        Q->nInfo.setCId(-1);
                        return;
                    }
                }
                Q->nInfo.setCId(id);
            }
        }

        inline intT SwitchMode(nodeT *Q, intT id){
            if(isTopDown(id)) return id;
            intT cid = Q->nInfo.getCId();
            if(cid == -1) return -1;
            return uf->find(cid); 
        }

        inline bool Par(nodeT *Q){return Q->n > 2000;}

        inline bool Stop(nodeT *Q, intT id){
            intT cid = Q->nInfo.getCId();
            return cid != -1 && cid == id;
        }

    };

    template<class nodeT, class F, class E>
    void singletree(nodeT *Q, F *f, E id){
        id = f->SwitchMode(Q, id);
        if(f->Stop(Q, id)) return;
        if(Q->isLeaf()){
            f->BaseCase(Q, id);
        }else{
            f->TopDownNode(Q, id);
            if(f->Par(Q)){
                cilk_spawn singletree<nodeT, F, E>(Q->left, f, id);
                singletree<nodeT, F, E>(Q->right, f, id);
                cilk_sync;
            }else{
                singletree<nodeT, F, E>(Q->left, f, id);
                singletree<nodeT, F, E>(Q->right, f, id);
            }
            f->BottomUpNode(Q, id);
        }
    }

    template<int dim, class nodeT>
    void markTree(nodeT *Q, UnionFind::ParUF<intT> *uf, intT cid = -1){
        typedef MarkClusterId<dim, nodeT> M;
        M *marker = new M(uf);
        singletree<nodeT, M, intT>(Q, marker, cid);
    }


    template<int dim, class nodeT>
    void printTreeCId(nodeT *Q){
        typedef UTIL::PrintClusterId<dim, nodeT> M;
        M *printer = new M();
        singletree<nodeT, M, intT>(Q, printer, -1);
    }
}

#endif