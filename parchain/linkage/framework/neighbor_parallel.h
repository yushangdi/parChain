#ifndef FINDNNP_H
#define FINDNNP_H

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include "shared.h"
#include "kdTree2.h"
#include "ndHash.h"
#include "bruteforce.h"
#include "neighbor.h"

using namespace std;

namespace FINDNNP {

#define dualtree_spawn_macro(_condA, _condB, _Q1, _R1, _Q2, _R2){ \
    if(_condA && _condB){ \
       {cilk_spawn dualtree<nodeT, F>(_Q1, _R1, f, false); \
        dualtree<nodeT, F>(_Q2, _R2, f, false); \
        cilk_sync; }\
    }else if(_condA){ \
        dualtree<nodeT, F>(_Q1, _R1, f, false); \
    }else{ \
        dualtree<nodeT, F>(_Q2, _R2, f, false); \
    } }

    template<class nodeT, class F>
    void dualtree(nodeT *Q, nodeT *R, F *f, bool check = true){
        if(Q->size() + R->size() < 4000){
            FINDNN::dualtree_serial<nodeT, F>(Q, R, f, check);
            return;
        }

        if(f->Score(Q,R,check)) return;

        if(f->isLeaf(Q) && f->isLeaf(R)){ // not parallel bc leaf size should be small
            // for(intT i=0; i<Q->n; ++i){
            //     for(intT j=0;j<R->n; ++j){
            //         f->BaseCase(Q,R,i,j);
            //     }
            // }
            f->BasePost(Q,R);
        }else if (f->isLeaf(Q)){
           double dLeft = f->NodeDistForOrder(Q, R->left);
	       double dRight = f->NodeDistForOrder(Q, R->right);
            bool condA = !f->Score(dRight, Q, R->right);
            bool condB = !f->Score(dLeft, Q, R->left);
           if(f->SpawnOrder(dLeft , dRight)){
               dualtree_spawn_macro(condA, condB, Q, R->right, Q, R->left);
           }else{
               dualtree_spawn_macro(condB, condA, Q, R->left, Q, R->right);
           }
           f->QLPost(Q,R);
        }else if(f->isLeaf(R)){
           double dLeft = f->NodeDistForOrder(Q->left, R);
	       double dRight = f->NodeDistForOrder(Q->right, R);
            bool condA = !f->Score(dRight,Q->right, R);
            bool condB = !f->Score(dLeft, Q->left, R);
           if(f->SpawnOrder(dLeft , dRight)){
               dualtree_spawn_macro(condA, condB, Q->right, R, Q->left, R);
           }else{
              dualtree_spawn_macro(condB, condA, Q->left, R, Q->right, R);
           }
            f->RLPost(Q,R);
        }else{
            pair<double, pair<nodeT *, nodeT *>> callOrder[4];
            callOrder[0] = make_pair(f->NodeDistForOrder(Q->left, R->left), make_pair(Q->left, R->left));
            callOrder[1] = make_pair(f->NodeDistForOrder(Q->right, R->left), make_pair(Q->right, R->left));
            callOrder[2] = make_pair(f->NodeDistForOrder(Q->left, R->right), make_pair(Q->left, R->right));
            callOrder[3] = make_pair(f->NodeDistForOrder(Q->right, R->right), make_pair(Q->right, R->right));
            sort(callOrder, callOrder + 4);
            // parallel_for_1 (int cc = 0; cc < 4; ++ cc) {//TODO: expand to spawn
            //     int c = f->SpawnOrder(cc);
            //     nodeT *QQ = callOrder[c].second.first;
            //     nodeT *RR = callOrder[c].second.second;
            //     if(!f->Score(callOrder[c].first, QQ, RR)) dualtree<nodeT, F>(QQ, RR, f, false);
            // }
                cilk_spawn dualtree<nodeT, F>(callOrder[f->SpawnOrder(0)].second.first, callOrder[f->SpawnOrder(0)].second.second, f, true);
                cilk_spawn dualtree<nodeT, F>(callOrder[f->SpawnOrder(1)].second.first, callOrder[f->SpawnOrder(1)].second.second, f, true);
                cilk_spawn dualtree<nodeT, F>(callOrder[f->SpawnOrder(2)].second.first, callOrder[f->SpawnOrder(2)].second.second, f, true);
                dualtree<nodeT, F>(callOrder[f->SpawnOrder(3)].second.first, callOrder[f->SpawnOrder(3)].second.second, f, true);
            cilk_sync;
            f->Post(Q,R);
        }
    }


    template<class nodeT, class F, class E>
    E dualtree2(nodeT *Q, nodeT *R, F *f, E &result, bool check = true){
        if(Q->size() + R->size() < 200){
            result = FINDNN::dualtree_serial2<nodeT, F, E>(Q, R, f, check);
            return result;
        }
        // E result;
        bool stop;
        tie(result, stop) = f->Score(Q,R,check);
        if(stop) return result;

        if(f->isLeaf(Q) && f->isLeaf(R)){
            result = f->BasePost(Q,R);
        }else if (f->isLeaf(Q)){
            E result1, result2;
            cilk_spawn dualtree2<nodeT, F, E>(Q, R->left, f, result1, false);
            result2 = dualtree2<nodeT, F, E>(Q, R->right, f, result2, false);
            cilk_sync;
            result = f->QLPost(Q,R, result1, result2);
        }else if(f->isLeaf(R)){
            E result1, result2;
            cilk_spawn dualtree2<nodeT, F, E>(Q->left, R, f, result1, false);
            result2 = dualtree2<nodeT, F, E>(Q->right, R, f, result2, false);
            cilk_sync;
            result = f->RLPost(Q,R, result1, result2);
        }else{
            E result1, result2, result3, result4;
            cilk_spawn dualtree2<nodeT, F, E>(Q->left, R->left, f, result1, false);
            cilk_spawn dualtree2<nodeT, F, E>(Q->left, R->right, f, result2, false);
            cilk_spawn dualtree2<nodeT, F, E>(Q->right, R->right, f, result3, false);
            result4 = dualtree2<nodeT, F, E>(Q->right, R->left, f, result4, false);
	        cilk_sync;
            result = f->Post(Q,R, result1, result2, result3, result4);
        }
        return result;
    }

}
#endif
