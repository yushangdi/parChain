#ifndef FINDNN_NODE_H
#define FINDNN_NODE_H

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include "shared.h"
#include "kdTree2.h"
using namespace std;

namespace FINDNN {

    struct dummyNodeInfo{};

    struct CLinkNodeInfo{
        typedef intT infoT; 
        intT cId;
        double ub;
        intT round = 0;
        intT idx = -1; // node idx

        CLinkNodeInfo(){
            cId = -1;
            ub = numeric_limits<double>::max(); // used for allptsnn
        }

        inline double getUB(){return ub;}
        inline void updateUB(double tmp){
            utils::writeMin(&ub, tmp);
        }
        inline void resetUB(){
            ub = numeric_limits<double>::max();
        }
        inline intT getCId(){return cId;}
        inline void setCId(intT id){cId = id;}
        inline intT getRound(){return round;}
        inline void setRound(intT r){round = r;}
        inline intT getIdx(){return idx;}
        inline void setIdx(intT r){idx = r;}
        inline infoT initInfoVal() {return -1;}
        inline intT getMinN(){return 1;}

    };

    struct HLinkNodeInfo: public CLinkNodeInfo {
        double max_core_dist = -1;
        HLinkNodeInfo(){
            cId = -1;
            ub = numeric_limits<double>::max(); // used for allptsnn
        }
    };


    struct WLinkNodeInfo{//TODO: remove cid
        typedef tuple<intT, intT> infoT; // <cid, min_n>
        intT min_n = 1;
        intT cId;
        double ub;
        // intT idx = -1; // node idx

        WLinkNodeInfo(){
            cId = -1;
            ub = numeric_limits<double>::max(); // used for allptsnn
        }

        inline infoT initInfoVal() {return make_tuple(-1, numeric_limits<intT>::max());}

        inline void setInfo(infoT info){
            tie(cId, min_n) = info;
        }
        inline infoT getInfo(){
            return infoT(cId, min_n);
        }
        inline intT getMinN(){
            return min_n;
        }
        inline void setMinN(intT nn){
            min_n = nn;
        }
        inline double getUB(){return ub;}
        inline void updateUB(double tmp){
            utils::writeMin(&ub, tmp);
        }
        inline intT getCId(){return cId;}
        inline void setCId(intT id){cId = id;}
        // inline intT getIdx(){return idx;}
        // inline void setIdx(intT r){idx = r;}
    };

    struct WLinkNodeInfo2{//TODO: remove cid
        typedef intT infoT; // min_n
        intT min_n = 1;
        // intT cId;
        double ub;
        // intT idx = -1; // node idx

        WLinkNodeInfo2(){
            // cId = -1;
            ub = numeric_limits<double>::max(); // used for allptsnn
        }

        inline infoT initInfoVal() {return numeric_limits<intT>::max();}

        inline void setInfo(infoT info){
            min_n = info;
        }
        inline infoT getInfo(){
            return min_n;
        }
        inline intT getMinN(){
            return min_n;
        }
        inline void setMinN(intT nn){
            min_n = nn;
        }
        inline double getUB(){return ub;}
        inline void updateUB(double tmp){
            utils::writeMin(&ub, tmp);
        }
        inline intT getCId(){return -1;}
        inline void setCId(intT id){;}
    };

    // linkedlist nodes
    template<int dim, class objT>
    struct LinkNodeInfo2{
        typedef point<dim> pointT;
        intT cId;
        intT round = 0;
        intT idx; // node idx
        intT n = 1;
        LDS::node<objT> *items;
        LinkNodeInfo2<dim, objT> *left = nullptr;
        LinkNodeInfo2<dim, objT> *right = nullptr;
        LDS::node<objT> *tail;
        pointT pMin, pMax;
        double height = 0;

        LinkNodeInfo2(intT  t_cid, intT t_round, intT t_idx, LinkNodeInfo2<dim, objT> *t_left, LinkNodeInfo2<dim, objT> *t_right, double _height):
            cId(t_cid),
            round(t_round),
            idx(t_idx),
            items(t_left->items),
            left(t_left),
            right(t_right),
            height(_height){
                n = t_left->n + t_right->n;
                tail = t_right->tail;
                left->tail->next = t_right->items;
                pMin = pointT(left->pMin.x);
                pMax = pointT(left->pMax.x);
                pMin.minCoords(right->pMin);
                pMax.maxCoords(right->pMax);
            }

        LinkNodeInfo2(intT  t_cid, LDS::node<objT> *t_items):
            cId(t_cid),
            idx(t_cid),
            items(t_items),
            tail(t_items){
                pMin = t_items->elt->p;
                pMax = t_items->elt->p;
            }

        inline bool isLeaf(){return n == 1;}
        inline intT size(){return n;}
        inline double getHeight(){return height;}
        inline intT getRound(){return round;}
        inline intT getIdx(){return idx;}
    };

    // consecutive points
    template<int dim, class objT>
    struct LinkNodeInfo1{
        intT cId;
        intT round = 0;
        intT idx; // node idx
        intT n = 1;
        LinkNodeInfo1<dim, objT> *left = nullptr;
        LinkNodeInfo1<dim, objT> *right = nullptr;
        point<dim> pMin, pMax;
        intT offset = -1; 
        double height = 0;
        point<dim> center;

        LinkNodeInfo1(intT  t_cid, intT t_round, intT t_idx, LinkNodeInfo1<dim, objT> *t_left, LinkNodeInfo1<dim, objT> *t_right, double _height):
            cId(t_cid),
            round(t_round),
            idx(t_idx),
            left(t_left),
            right(t_right),
            height(_height){
                intT nl = t_left->n;
                intT nr = t_right->n;
                n = nl+nr;
                pMin = left->pMin;//pointT(left->pMin.x);
                pMax = left->pMax;//pointT(left->pMax.x);
                pMin.minCoords(right->pMin);
                pMax.maxCoords(right->pMax);
                if(height ==0){
                    center = t_left->center; //important to avoid numerical error
                }else{
                    for (int i=0; i<dim; ++i) {
                        center.updateX(i, (nl * t_left->center[i] + nr * t_right->center[i]) / n);
                    }
                }
            }
        
        LinkNodeInfo1(intT  t_cid, iPoint<dim> p):
            cId(t_cid),
            idx(t_cid),
            offset(t_cid){
                pMin = p.p;
                pMax = p.p;
                center = p.p;
            }

        inline bool isLeaf(){return n == 1;}
        inline intT getOffset(){return offset;}
        inline void setOffset(intT i){offset = i;}
        inline double getHeight(){return height;}
        inline intT getRound(){return round;}
        inline intT getIdx(){return idx;}
        inline intT size(){return n;}
        inline double dist(LinkNodeInfo1<dim, objT> *node){return center.pointDistSq(node->center);}
    };

    // consecutive points
    template<int dim, class objT>
    struct LinkNodeInfo3{
        intT cId;
        intT round = 0;
        intT idx; // node idx
        intT n = 1;
        LinkNodeInfo3<dim, objT> *left = nullptr;
        LinkNodeInfo3<dim, objT> *right = nullptr;
        point<dim> center;
        double height = 0;
        // point<dim> pMin, pMax;
        // intT offset = -1; 

        LinkNodeInfo3(intT  t_cid, intT t_round, intT t_idx, LinkNodeInfo3<dim, objT> *t_left, LinkNodeInfo3<dim, objT> *t_right, double _height):
            cId(t_cid),
            round(t_round),
            idx(t_idx),
            left(t_left),
            right(t_right),
            height(_height){
                intT nl = t_left->n;
                intT nr = t_right->n;

                n = t_left->n + t_right->n;
                if(height ==0){
                    center = t_left->center; //important to avoid numerical error
                }else{
                    for (int i=0; i<dim; ++i) {
                        center.updateX(i, (nl * t_left->center[i] + nr * t_right->center[i]) / n);
                    }
                }
        }
        
        LinkNodeInfo3(intT  t_cid, iPoint<dim> p):
            cId(t_cid),
            idx(t_cid){
                center = p.p;
            }

        inline bool isLeaf(){return n == 1;}
        inline double dist(LinkNodeInfo3<dim, objT> *node){return center.pointDistSq(node->center);}
        inline double getHeight(){return height;}
        inline intT getRound(){return round;}
        inline intT getIdx(){return idx;}
        inline intT size(){return n;}
    };

}//end  FINDNN
#endif