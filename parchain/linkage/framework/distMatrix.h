#ifndef METHOD_CHAIN_TREE_DM_H
#define METHOD_CHAIN_TREE_DM_H

#include <limits>
#include <math.h>
#include <string.h>

#include <string>
#include <vector>
#include <stdio.h>
#include <algorithm>


#include "geometry.h"
#include "shared.h"
#include "dist.h"

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

namespace MatrixDistanceComputer{

template<int dim, class pointTT, class nodeTT>
struct distMatrixAbstract: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;
  pointT *PP;

  pointT *centers; //ith location is the center of cluster activeCluster[i], need to be contiguous when building tree
  bool use_range;


  distMatrixAbstract(pointT *t_PP, intT t_n, bool t_use_range): PP(t_PP){
    use_range = t_use_range;

    if(use_range){
    centers = newA(pointT, t_n);
    parallel_for(intT i=0; i<t_n; ++i) {centers[i] = PP[i];}
    }

  }

  inline void initNodes(nodeT *nodes, intT n){
    parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
  }

  template<class F>
  inline void update(intT round, F *finder){

      if(use_range){
    //make points array
    intT  C = finder->C;

    parallel_for(intT i = 0; i < C; ++i){
      intT cid = finder->activeClusters[i];
      nodeT *clusterNode = finder->getNode(cid);
      centers[i] = pointT(clusterNode->center, cid);
    }

    // build kdtree
    finder->kdtree->kdTreeRebuild(centers, C);
      }

  }

  inline double updateDistO(double d1, double d2, double nql, double nqr, double nr, double dij){
      return 0;
  }

  inline double updateDistN(double d1, double d2, double d3, double d4, 
                       double nql, double nqr, double nrl, double nrr,
                       double dij, double dklr){
    return 0;
  }

  ~distMatrixAbstract(){
    if(use_range)free(centers);
  }

};

template<int dim, class pointTT, class nodeTT, class nodeInfo>
struct distAverage4: public distMatrixAbstract<dim, pointTT, nodeTT> {
  typedef pointTT pointT;
  typedef nodeTT nodeT;

  bool id_only = true;
  LDS::Method method = LDS::AVGSQ;
  inline static void printName(){
    cout << "dist Avg 4, use centers and variances, norm = 2" << endl;
  }
    
  distAverage4(pointT *t_PP, intT t_n, bool t_use_range): distMatrixAbstract<dim, pointTT, nodeTT>(t_PP, t_n, t_use_range){
  }

  inline double updateDistO(double d1, double d2, double nql, double nqr, double nr, double dij){
    double n1 = (double)nql * (double)nr;
    double n2 = (double)nqr * (double)nr;
    double alln = n1 + n2 ;
    d1 = n1 / alln * d1;
    d2 = n2 / alln * d2;
    return d1 + d2;
  }

  inline double updateDistN(double d1, double d2, double d3, double d4, 
                       double nql, double nqr, double nrl, double nrr,
                       double dij, double dklr){
    double n1 = (double)nql * (double)nrl;
    double n2 = (double)nql * (double)nrr;
    double n3 = (double)nqr * (double)nrl;
    double n4 = (double)nqr * (double)nrr;
    double alln = n1 + n2 + n3 + n4;
    d1 = n1 / alln * d1;
    d2 = n2 / alln * d2;
    d3 = n3 / alln * d3;
    d4 = n4 / alln * d4;
    return d1  + d2  + d3 + d4;
  }

  template<class F>
  inline void postProcess(F *finder){
  }

  ~distAverage4(){
  }

};

// group clustered points at the end of each round
// compute distance by traverse clustered points
// same as distAverage1, except rebuild kdtree at the end
template<int dim, class pointTT, class nodeTT>
struct distAverage5: public distMatrixAbstract<dim, pointTT, nodeTT> {
  typedef pointTT pointT;
  typedef nodeTT nodeT;

  LDS::Method method = LDS::AVG;

  inline static void printName(){
    cout << "distAverage5 point array, norm = 1, rebuild tree" << endl;
  }

  distAverage5(pointT *t_PP, intT t_n, bool t_use_range): distMatrixAbstract<dim, pointTT, nodeTT>(t_PP, t_n, t_use_range){
  }

  
  inline double updateDistO(double d1, double d2, double nql, double nqr, double nr, double dij){
    double n1 = (double)nql * (double)nr;
    double n2 = (double)nqr * (double)nr;
    double alln = n1 + n2 ;
    d1 = n1 / alln * d1;
    d2 = n2 / alln * d2;
    return d1 + d2;
  }

  inline double updateDistN(double d1, double d2, double d3, double d4, 
                       double nql, double nqr, double nrl, double nrr,
                       double dij, double dklr){
    double n1 = (double)nql * (double)nrl;
    double n2 = (double)nql * (double)nrr;
    double n3 = (double)nqr * (double)nrl;
    double n4 = (double)nqr * (double)nrr;
    double alln = n1 + n2 + n3 + n4;
    d1 = n1 / alln * d1;
    d2 = n2 / alln * d2;
    d3 = n3 / alln * d3;
    d4 = n4 / alln * d4;
    return d1  + d2  + d3 + d4;
  }

  ~distAverage5(){
  }

};

//////////////////////////////////// Complete Linkage ///////////////////////////////

// same as distComplete1, uses array instead of cache
template<int dim, class pointTT, class nodeTT, class nodeInfo>
struct distComplete3: public distMatrixAbstract<dim, pointTT, nodeTT> {
  typedef pointTT pointT;
  typedef nodeTT nodeT;
  typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
  typedef kdTree<dim, pointT, nodeInfo> kdtreeT;
  typedef FINDNN::NNcomplete1<dim, kdnodeT> FComp;

  typedef tuple<intT, intT, long> clusterCacheT;//(round, count, cid), use long to pack to 2^i bytes
  // typedef pair<pair<intT, intT>, pair<intT, intT>> clusterCacheT; 

  LDS::Method method = LDS::COMP;

  inline static void printName(){
    cout << "distComplete3 point* array, kdtree, array for counting" << endl;
  }


  distComplete3(pointT *t_PP, intT t_n, bool t_use_range): distMatrixAbstract<dim, pointTT, nodeTT>(t_PP, t_n, t_use_range){
  }

  inline double updateDistO(double d1, double d2, double nql, double nqr, double nr, double dij){
    return max(d1,d2);
  }

  inline double updateDistN(double d1, double d2, double d3, double d4, 
                       double nql, double nqr, double nrl, double nrr,
                       double dij, double dklr){
    return max(max(max(d1,d2), d3), d4);
  }


  ~distComplete3(){
  }

};



//////////////////////////////////// Ward's Linkage ///////////////////////////////

template<int dim, class pointTT, class nodeTT, class nodeInfo, class M>
struct distWard1: public distMatrixAbstract<dim, pointTT, nodeTT> {
  typedef pointTT pointT;
  typedef nodeTT nodeT;
  typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
  typedef kdTree<dim, pointT, nodeInfo> kdtreeT;

  LDS::Method method = LDS::WARD;
  pointT *PP;
  M marker;
  intT *sizes; //ith location is the size of cluster i

  inline static void printName(){
    cout << "dist Ward 1" << endl;
  }
    
  distWard1(pointT *t_PP, intT t_n, bool t_use_range): distMatrixAbstract<dim, pointTT, nodeTT>(t_PP, t_n, t_use_range){
      if(t_use_range){
        sizes = newA(intT, t_n);
        parallel_for(intT i=0; i<t_n; ++i) {sizes[i] = 1;}
        marker = M(sizes);
      }
    
  }
  
  template<class F>
  inline void update(intT round, F *finder){
      if(this->use_range){
    //make points array
    intT  C = finder->C;
    
    parallel_for(intT i = 0; i < C; ++i){
      intT cid = finder->activeClusters[i];
      this->centers[i] = pointT(finder->getNode(cid)->center, cid);
    }

    // build kdtree
    finder->kdtree->kdTreeRebuild(this->centers, C);

    //mark min sizes
    FINDNN::singletree<kdnodeT, M, typename M::infoT>(finder->kdtree->root, &marker, marker.initVal);
      }

  }

  inline double updateDistO(double dik, double djk, double ni, double nj, double nk, double dij){
    double ntotal = ni + nj + nk;
    double d = sqrt( ( ((ni + nk)  * dik * dik) + ((nj + nk) * djk * djk) - (nk * dij * dij) )/ ntotal );
    return d;
  }

  inline double updateDistN(double dikl, double dikr, double djkl, double djkr, 
                       double ni, double nj, double nkl, double nkr,
                       double dij, double dklr ){
    double dik = updateDistO(dikl, dikr, nkl, nkr, ni, dklr);
    double djk = updateDistO(djkl, djkr, nkl, nkr, nj, dklr);
    return updateDistO(dik, djk, ni, nj, (nkl + nkr), dij);
  }

  ~distWard1(){
      if(this->use_range) free(sizes);
  }

};
}

#endif