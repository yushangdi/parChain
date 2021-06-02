#ifndef CLINK_BRUTEFORCENN_H
#define CLINK_BRUTEFORCENN_H

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include "shared.h"
#include "kdTree2.h"
#include "unionfind.h"

// 10000 similar runtime
#define LINKAGE_AVE_BRUTE_THRESH 100000

#define BruteStrategy4

// Bruteforce implementations
namespace FINDNN {

#if defined(BruteStrategy4)
    //find the average　distance between P1 and P2, return random pair of points
    //TODO: tune
    template<class pointT, class ARR>
    inline pair<pair<intT, intT>, double> bruteForceAverage(ARR P1, ARR P2, intT n1, intT n2, bool div = true){
        pair<intT, intT> result = make_pair(LDS::getPointId(P1[0]),LDS::getPointId(P2[0]));
        double total_d = 0;
        long total_n = (long)n1 * (long)n2;
        // intT eltsPerCacheLine = 128 /sizeof(double);
        if(total_n < LINKAGE_AVE_BRUTE_THRESH){ //threshold?
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

  // template <class OT, class intT, class F, class G>
  // OT reduce2(intT s, intT e, F f, G g) {
  //   timer t1; t1.start();
  //   intT nBlocks1 =  getWorkers() * 32;
  //   intT BLOCKSIZE1 = max((intT)1,(e-s)/nBlocks1);
  //   intT l = nblocks(e-s, BLOCKSIZE1);
  //   if (l <= 1) return reduceSerial<OT>(s, e, f , g);
  //   OT *Sums = newA(OT,l);
  //   // blocked_for (i, s, e, BLOCKSIZE1,
  //   //  Sums[i] = reduceSerial<OT>(s, e, f, g););
  //   cout << "reduce part 0: " << t1.next() << endl;
  //   cout << l << endl;
  //   parallel_for_1 (intT _i = 0; _i < l; _i++) {   
  //     intT ss = s + _i * BLOCKSIZE1;      
  //     intT ee = min(ss + BLOCKSIZE1, e);      
  //     Sums[_i] = reduceSerial<OT>(ss, ee, f, g);        
  //   }     
  //   cout << "reduce part 1: " << t1.next() << endl;
  //   OT r = reduce<OT>((intT) 0, l, f, getA<OT,intT>(Sums));
  //   free(Sums);
  //   return r;
  // }


#elif defined(BruteStrategy5)
  template <class ARR, class pointT>
  struct getABruteAve {
    ARR P1;
    ARR P2;
    // intT n1;
    long n2;
    // long shift;
    // long mask;
    getABruteAve(ARR PP1, ARR PP2, long nn2) : P1(PP1), P2(PP2), n2(nn2) {
        // shift = log2(n2);
        // mask = n2-1;
    }
    double operator() (long k) {
        long i = k/n2;
        long j = k - n2*i;
        // long i = k >> shift;
        // long j = k & mask;
        return LDS::pointDist<pointT>(P1[i], P2[j]);
    }
  };


    //find the average　distance between P1 and P2, return random pair of points
    template<class pointT, class ARR>
    inline pair<pair<intT, intT>, double> bruteForceAverage(ARR P1, ARR P2, intT n1, intT n2, bool div = true){
        pair<intT, intT> result = make_pair(LDS::getPointId(P1[0]),LDS::getPointId(P2[0]));
        double total_d = 0;
        long total_n = (long)n1 * (long)n2;
        // intT eltsPerCacheLine = 128 /sizeof(double);
        if(total_n < LINKAGE_AVE_BRUTE_THRESH){ //threshold?
            for(intT i=0; i< n1; ++i){
                for(intT j=0;j<n2; ++j){
                    total_d +=  LDS::pointDist<pointT>(P1[i], P2[j]);
                }
            }
        }else{

            total_d = sequence::reduce<double, long, utils::addF<double>, getABruteAve<ARR, pointT>>(
                (long)0, total_n, utils::addF<double>(),  getABruteAve<ARR, pointT>(P1, P2, (long)n2));

        }

        if(div) return make_pair(result, total_d/total_n);
        return make_pair(result, total_d);
    }


#elif defined(BruteStrategy3)
    template<class pointT, class ARR>
    inline pair<pair<intT, intT>, double> bruteForceAverage(ARR P1, ARR P2, intT n1, intT n2, bool div = true){
        pair<intT, intT> result = make_pair(LDS::getPointId(P1[0]),LDS::getPointId(P2[0]));
        double total_d = 0;
        long total_n = (long)n1 * (long)n2;
        intT eltsPerCacheLine = 128 /sizeof(double);
        if(total_n < LINKAGE_AVE_BRUTE_THRESH){ //threshold?
            for(intT i=0; i< n1; ++i){
                for(intT j=0;j<n2; ++j){
                    total_d += LDS::pointDist<pointT>(P1[i], P2[j]);
                }
            }
        }else{
            // long BLOCKSIZE = 2048;
            long nBlocks =  getWorkers();
            // long BLOCKSIZE = total_n/nBlocks;
            // if(total_n < nBlocks*BLOCKSIZE) nBlocks += 1;
            double *arr = newA(double, nBlocks * eltsPerCacheLine);
            for(intT i=0; i< nBlocks; ++i){
                arr[i*eltsPerCacheLine] = 0;
            }
            parallel_for(intT i=0; i< n1; ++i){
                parallel_for(intT j=0; j< n2; ++j){
                    arr[getWorkerId() * eltsPerCacheLine] += LDS::pointDist<pointT>(P1[i], P2[j]);
                }
            }
            for(intT i = 0; i < nBlocks; ++i){
                total_d += arr[i * eltsPerCacheLine];
            }
            // total_d = sequence::plusReduce(arr, nBlocks);
            free(arr);
        }

        if(div) return make_pair(result, total_d/total_n);
        return make_pair(result, total_d);
    }
#endif

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
            // else if(d==maxd && j<result.second){
            //     result = make_pair(i,j);
            //     maxd = d;
            // }
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
            // else if(d==maxd && P2[j]->idx()<result.second){
            //     result = make_pair(P1[i]->idx(),P2[j]->idx());
            //     maxd = d;
            // }
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
