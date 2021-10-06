#ifndef CLINK_BRUTEFORCENN_H
#define CLINK_BRUTEFORCENN_H

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include "shared.h"
#include "kdTree2.h"
#include "unionfind.h"

#define LINKAGE_AVE_BRUTE_THRESH 100000

// Bruteforce implementations
namespace FINDNN {

    //find the average　distance between P1 and P2, return random pair of points
    template<class pointT, class ARR>
    inline pair<pair<intT, intT>, double> bruteForceAverage(ARR P1, ARR P2, intT n1, intT n2, bool div = true){
        pair<intT, intT> result = make_pair(LDS::getPointId(P1[0]),LDS::getPointId(P2[0]));
        double total_d = 0;
        long total_n = (long)n1 * (long)n2;
        // intT eltsPerCacheLine = 128 /sizeof(double);
        if(total_n < LINKAGE_AVE_BRUTE_THRESH){ 
            for(intT i=0; i< n1; ++i){
                for(intT j=0;j<n2; ++j){
                    total_d +=  LDS::pointDist<pointT>(P1[i], P2[j]);
                }
            }
        }else{
            if(n1 < n2){
                swap(n1, n2);
                swap(P1, P2);
            }
            intT nBlocks1 =  getWorkers() * 8;

            intT BLOCKSIZE1 = max(1,n1/nBlocks1);

            if(n1 > nBlocks1*BLOCKSIZE1) BLOCKSIZE1 += 1;

            parallel_for_1(intT ii=0; ii< nBlocks1; ++ii){
                double thread_sum = 0;
                for (intT i=ii*BLOCKSIZE1; i< min((ii+1)*BLOCKSIZE1, n1); ++i){
                    for (intT j=0; j< n2; ++j){
                        thread_sum +=  LDS::pointDist<pointT>(P1[i], P2[j]);
                    }
                }
                utils::writeAdd(&total_d, thread_sum);
            }

        }

        if(div) return make_pair(result, total_d/total_n);
        return make_pair(result, total_d);
    }

    template<int dim, class pointT, class nodeInfo>
    inline pair<pair<intT, intT>, double> bruteForceAverage(struct kdTree<dim, pointT, nodeInfo> *t1, struct kdTree<dim, pointT, nodeInfo> *t2, bool div = true){
        pointT** P1 = t1->items;
        pointT** P2 = t2->items;
        intT n1 = t1->getN();
        intT n2 = t2->getN();
        return bruteForceAverage<pointT, pointT**>(P1, P2, n1, n2, div);
    }

    template<int dim, class pointT, class nodeInfo>
    inline pair<pair<intT, intT>, double> bruteForceAverage(struct kdNode<dim, pointT, nodeInfo> *t1, struct kdNode<dim, pointT, nodeInfo> *t2, bool div = true){
        pointT** P1 = t1->items;
        pointT** P2 = t2->items;
        intT n1 = t1->n;
        intT n2 = t2->n;
        return bruteForceAverage<pointT, pointT**>(P1, P2, n1, n2, div);
    }

    template<class pointT, class nodeT>
    inline pair<pair<intT, intT>, double> bruteForceAverage(nodeT *inode,  nodeT *jnode, pointT **clusteredPts){
        pointT **P1 = clusteredPts + inode->getOffset();
        pointT **P2 = clusteredPts + jnode->getOffset();
        intT n1 = inode->n;
        intT n2 = jnode->n;
        return FINDNN::bruteForceAverage<pointT, pointT**>(P1, P2, n1, n2, true);
    }

    template<class pointT, class nodeT>
    inline pair<pair<intT, intT>, double> bruteForceAverage(nodeT *inode,  nodeT *jnode, pointT *clusteredPts){
        pointT *P1 = clusteredPts + inode->getOffset();
        pointT *P2 = clusteredPts + jnode->getOffset();
        intT n1 = inode->n;
        intT n2 = jnode->n;
        return FINDNN::bruteForceAverage<pointT, pointT*>(P1, P2, n1, n2, true);
    }

    template<int dim>
    inline pair<pair<intT, intT>, double> bruteForceFarthest(point<dim>* P, intT start, intT n){
        pair<intT, intT> result = make_pair(-1,-1);
        double maxd = 0;
        for(intT i=0; i< start; ++i){
            for(intT j=start;j<n; ++j){
            double d = P[i].pointDist(P[j]);
            if(d > maxd){
                result = make_pair(i,j);
                maxd = d;
            }
            }
        }
        return make_pair(result, maxd);
    }

    //find the farthest pair of points in t1 and t2
    template<class pointT>
    inline pair<pair<intT, intT>, double> bruteForceFarthest(pointT** P1, pointT** P2, intT n1, intT n2){
        pair<intT, intT> result = make_pair(-1,-1);
        double maxd = 0;
        for(intT i=0; i< n1; ++i){
            for(intT j=0;j<n2; ++j){
            double d = P1[i]->pointDist(P2[j]);
            if(d > maxd){
                result = make_pair(P1[i]->idx(),P2[j]->idx());
                maxd = d;
            }
            }
        }
        return make_pair(result, maxd);
    }

    template<int dim, class pointT, class nodeInfo>
    inline pair<pair<intT, intT>, double> bruteForceFarthest(struct kdTree<dim, pointT, nodeInfo> *t1, struct kdTree<dim, pointT, nodeInfo> *t2){
        pointT** P1 = t1->items;
        pointT** P2 = t2->items;
        intT n1 = t1->getN();
        intT n2 = t2->getN();
        return bruteForceFarthest(P1, P2, n1, n2);
    }

    template<class pointT>
    inline pair<pair<intT, intT>, double> bruteForceNearest(pointT** P1, pointT** P2, intT n1, intT n2){
        pair<intT, intT> result;
        double mind = numeric_limits<double>::max();
        for(intT i=0; i< n1; ++i){
            for(intT j=0;j<n2; ++j){
            double d = P1[i]->pointDist(P2[j]);
            if(d < mind){
                result = make_pair(P1[i]->idx(),P2[j]->idx());
                mind = d;
            }else if(d==mind && P2[j]->idx()<result.second){
                result = make_pair(P1[i]->idx(),P2[j]->idx());
                mind = d;
            }
            }
        }
        return make_pair(result, mind);
    }

    template<class pointT>
    inline pair<pair<intT, intT>, double> bruteForceNearest(pointT* P1, pointT* P2, intT n1, intT n2){
        pair<intT, intT> result;
        double mind = numeric_limits<double>::max();
        for(intT i=0; i< n1; ++i){
            for(intT j=0;j<n2; ++j){
            double d = P1[i].pointDist(P2[j]);
            if(d < mind){
                result = make_pair(P1[i].idx(),P2[j].idx());
                mind = d;
            }else if(d==mind && P2[j].idx()<result.second){
                result = make_pair(P1[i].idx(),P2[j].idx());
                mind = d;
            }
            }
        }
        return make_pair(result, mind);
    }



    //find the farthest pair of points in t1 and t2
    template<int dim, class pointT, class nodeInfo>
    inline pair<pair<intT, intT>, double> bruteForceNearest(struct kdTree<dim, pointT, nodeInfo> *t1, struct kdTree<dim, pointT, nodeInfo> *t2){
        pointT** P1 = t1->items;
        pointT** P2 = t2->items;
        intT n1 = t1->getN();
        intT n2 = t2->getN();

        return bruteForceNearest(P1, P2, n1, n2);
    }

    template<int dim, class pointT>
    inline bool itemInBox(pointT pMin1, pointT pMax1, pointT* item) {
    for(int i=0; i<dim; ++i) {
      if (pMax1[i]<item->x[i] || pMin1[i]>item->x[i]) return false;
    }
    return true;}

    //find all points within 2r side-length blocks centered at query
    // return cluster ids and number of points in each cluster
    template<int dim>
    inline intT* bruteForceNNCandidateBox(point<dim> query, point<dim>* P, intT n, double r, UnionFind::ParUF<intT> *uf){
        intT  *sizes = newA(intT, n);
        for(intT i=0; i< n; ++i){
            sizes[i] = 0;
        }

        point<dim> pMin1 = point<dim>();
        point<dim> pMax1 = point<dim>();
        double* center = query.x;
        for (int i=0; i<dim; ++i) {
        pMin1.updateX(i, center[i]-r);
        pMax1.updateX(i, center[i]+r);}

        for(intT i=0; i< n; ++i){
            if(itemInBox<dim, point<dim>>(pMin1, pMax1, P+i)){
                sizes[uf->find(i)]++;
            }
        }

        return sizes;
    }




}
#endif
