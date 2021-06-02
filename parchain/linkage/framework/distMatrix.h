#ifndef METHOD_CHAIN_TREE_DM_H
#define METHOD_CHAIN_TREE_DM_H

#include <limits>
#include <math.h>
#include <string.h>

#include <string>
#include <vector>
#include <stdio.h>
#include <algorithm>

// #include "sequence.h"
#include "geometry.h"
// #include "sampleSort.h"
// #include "gettime.h"
// #include "unionfind.h"
#include "shared.h"

// #include "../parallelHclust/fastcluster.cpp"
// #include "parseCommandLine.h"

template<int dim>
struct DM{

    intT n = 0;
    double* distmat = nullptr;
    long distmat_size = 0;

 DM(){}
 DM(point<dim>* P, intT _n, bool sqFlag = false){
    n = _n;

    timer t1;t1.start();
    // computation of condensed distance matrix
    distmat_size = (long)n*(n-1)/2;
    distmat = newA(double, distmat_size+1);//(double *)malloc(distmat_size * sizeof(double));//
    
    if(sqFlag){
        parallel_for (intT i=0; i<n; i++) {
            parallel_for (intT j=i+1; j<n; j++) {
                distmat[getInd(i, j)] = P[i].pointDistSq(P[j]);
            }
        }
    }else{
        parallel_for (intT i=0; i<n; i++) {
            parallel_for (intT j=i+1; j<n; j++) {
                distmat[getInd(i, j)] = P[i].pointDist(P[j]);
            }
        }
    }
    distmat[distmat_size] = 0;
    std::cout << "compute distance matrix " <<  t1.next() << std::endl;
 }

    inline long getInd(intT i, intT j){
        if(i == j) return distmat_size;
        // return (((2*distmat_size-3-(r_))*(r_))>>1)+(c_)-1;
        long r_ = static_cast<long>(i);
        long c_ = static_cast<long>(j);
        return (((2*n-r_-3) * r_) >> 1 )+ (c_)-1;
    }
    
 inline double get(intT r_, intT c_){
    if(r_ == c_) return 0;
    if(r_ > c_) swap(r_,c_);
    return( distmat[getInd(r_, c_)] );
 }

 inline void update(intT r_, intT c_, double dist){
    if(r_ == c_) return;
    if(r_ > c_) swap(r_,c_);
    distmat[getInd(r_, c_)] = dist;
 }

 ~DM(){
    free(distmat);
 }
    void printMatrix(){
    for (intT i=0; i<distmat_size; i++){
        cout << distmat[i] << endl;
    }
    for (intT i=0; i<n; i++){
        for (intT j=i+1; j<n; j++) {
            cout << i << " " << j << " " <<getInd(i,j) << " " << get(i,j) << endl;
        }
    } 
    }
};

//used for max index, so reverse comparason
struct neighborComparator{
    double eps;
    neighborComparator( double _eps):eps(_eps){}
    bool operator() (pair<double, intT> i, pair<double, intT> j) {
        double dist1 = i.first;
        double dist2 = j.first;
        if(abs(dist1 - dist2) <= eps) return i.second < j.second;
        return dist1 < dist2;
    }

};

template<class DM>
struct neighborDistGetter{
    intT cid;
    DM *matrix;
    intT *activeClusters;
    neighborDistGetter(intT _cid, DM *_matrix, intT *_activeClusters):cid(_cid), matrix(_matrix), activeClusters(_activeClusters){}
    pair<double, intT> operator() (intT i) {
        i = activeClusters[i];
        if(i==cid) return make_pair(numeric_limits<double>::max(), i); //used for max index, so reverse comparason
        return make_pair(matrix->get(i, cid), i);
    }

};

// template<int dim>
// int parallelHClust(point<dim>* P, intT n, commandLine params){
//     char* method = params.getOptionValue("-m");
//     int opt_method;
//     if (0 == strcmp(method, "single"))
//         opt_method = HCLUST_METHOD_SINGLE;
//     else if (0 == strcmp(method, "complete"))
//         opt_method = HCLUST_METHOD_COMPLETE;
//     else if (0 == strcmp(method, "average"))
//         opt_method = HCLUST_METHOD_AVERAGE;
//     else if (0 == strcmp(method, "median"))
//         opt_method = HCLUST_METHOD_MEDIAN;
//     else if (0 == strcmp(method, "ward"))
//         opt_method = HCLUST_METHOD_WARD;
//     else {
//         fputs("not a valid method", stderr);
//         return 1;
//     }


//     timer t1;t1.start();
//     intT npoints = n;
//     long distmat_size = ((long)n/2*(n-1));
//     double* distmat = new double[distmat_size];//(double *)malloc(distmat_size * sizeof(double));//
//     long k = 0;
//     parallel_for (intT i=0; i<n; i++) {
//         parallel_for (intT j=i+1; j<n; j++) {
//             distmat[k] = P[i].pointDist(P[j]);
//             k++;
//         }
//     }
//     std::cout << "compute distance matrix " <<  t1.next() << std::endl;
    
//     // clustering call
//     int* merge = new int[2*(npoints-1)];
//     double* height = new double[npoints-1];
//     hclust_fast(npoints, distmat, opt_method, merge, height);
//     std::cout << "hierarchy " <<  t1.next() << std::endl;
    
//     // int* labels = new int[npoints];
//     // cutree_k(npoints, merge, 2, labels);
//     // std::cout << "cut " <<  t1.next() << std::endl;
//     //cutree_cdist(npoints, merge, height, 0.5, labels);

//     double cost = 0;
//     for (intT i=0; i<npoints-1; i++) {
//       cost += height[i];
//     }
//     std::cout << std::setprecision(10) << "cost " <<  cost << std::endl;
    
//     // clean up
//     delete[] distmat;
//     delete[] merge;
//     delete[] height;
//     return 1;

// }

// template int parallelHClust<2>(point<2>*, intT, commandLine);
// template int parallelHClust<3>(point<3>*, intT, commandLine);
// template int parallelHClust<4>(point<4>*, intT, commandLine);
// template int parallelHClust<5>(point<5>*, intT, commandLine);
// template int parallelHClust<6>(point<6>*, intT, commandLine);
// template int parallelHClust<7>(point<7>*, intT, commandLine);
// template int parallelHClust<8>(point<8>*, intT, commandLine);
// template int parallelHClust<9>(point<9>*, intT, commandLine);
// template int parallelHClust<10>(point<10>*, intT, commandLine);
// template int parallelHClust<16>(point<16>*, intT, commandLine);
// template int parallelHClust<128>(point<128>*, intT, commandLine);

#endif