#ifndef LINKAGE_DIST_AVE_H
#define LINKAGE_DIST_AVE_H

using namespace std;
#include "neighbor.h"
#include "neighbor_parallel.h"

  struct distAbstract{

    bool id_only = false;
    bool nn_process = false;

    distAbstract(){}

    double getDistNaive(intT cid1, intT cid2, 
                        double lb = -1, double ub = numeric_limits<double>::max(), 
                        bool par = true){ 
      cout << "id only not supported" << endl; exit(1); return -1;
    }
    
    template<class kdnodeT, class Fs>
    inline void getRadius(intT cid, kdnodeT *root, Fs *fs){
      cout << "get Radius not supported. dist.h" << endl; exit(1);
    }

    template<class F>
    inline void postProcess(F *finder){
    }
  };



// group clustered points at the end of each round
// compute distance by traverse clustered points
template<int dim, class pointTT, class nodeTT>
struct distAverage1: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;

  // private:
  pointT *clusteredPts1;
  pointT *clusteredPts2;//clusteredPts point to one of the two
  pointT *PP;//used to initNodes only
  intT *clusterOffsets;  //same order as activeClusters in finder
  pointT *clusteredPts;

  // copy points from oldArray[oldOffset:oldOffset+n] to newArray[newOffset:newOffset+n]
  inline void copyPoints(pointT *oldArray, pointT *newArray, intT copyn, intT oldOffset, intT newOffset){
      parallel_for(intT j = 0; j < copyn; ++j){
          newArray[newOffset + j] = oldArray[oldOffset+j];
      }
    }
  
  inline static void printName(){
    cout << "distAverage1 point array, norm = 1" << endl;
  }

  distAverage1(pointT *t_PP, intT n, bool t_no_cache):PP(t_PP){
    clusteredPts1 = newA(pointT, n);
    clusteredPts2 = newA(pointT, n);
    clusteredPts = clusteredPts1;
    parallel_for(intT i=0; i<n; ++i) {
      clusteredPts1[i]=PP[i];
    }
    clusterOffsets = newA(intT, n);
    parallel_for(intT i=0; i<n; ++i) {clusterOffsets[i] = i;}
  }

    inline void initNodes(nodeT *nodes, intT n){
      parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
    }

    double getDistNaive(intT cid1, intT cid2, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
        return this->getDistNaive(cid1,cid2);
    }

    double getDistNaive(nodeT *inode,  nodeT *jnode, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    pair<pair<intT, intT>, double> result = FINDNN::bruteForceAverage(inode, jnode, clusteredPts);
    return result.second;
    }
  
    template<class kdnodeT2, class Fs>
    inline void getRadius(intT cid, kdnodeT2 *root, Fs *fs){
      this->getRadius(cid, root, fs);
    }



  template<class F>
  inline void update(intT round, F *finder){
    intT *activeClusters = finder->activeClusters;
    intT  C = finder->C;
    //  put points into clusters
    parallel_for(intT i = 0; i < C; ++i){
      clusterOffsets[i] = finder->getNode(activeClusters[i])->n;
    }
    sequence::prefixSum(clusterOffsets, 0, C);
    pointT *oldArray = clusteredPts;
    pointT *newArray = clusteredPts1;
    if(clusteredPts == clusteredPts1){
      newArray  = clusteredPts2;
    }

    //TODO: optimize for 2 clusters
    parallel_for(intT i = 0; i < C; ++i){
      intT cid  = activeClusters[i];
      nodeT *clusterNode = finder->getNode(cid);
      intT oldOffset = clusterNode->getOffset(); // must save this before setting new
      intT newOffset = clusterOffsets[i];
      clusterNode->setOffset(newOffset);
      if(clusterNode->round == round){//merged this round, copy left and right
        nodeT *clusterNodeL = clusterNode->left;
        oldOffset = clusterNodeL->getOffset();
        copyPoints(oldArray, newArray, clusterNodeL->n, oldOffset, newOffset);
        newOffset += clusterNodeL->n;

        nodeT *clusterNodeR = clusterNode->right;
        oldOffset = clusterNodeR->getOffset();
        copyPoints(oldArray, newArray, clusterNodeR->n, oldOffset, newOffset);

      }else{//not merged this round, just copy
        copyPoints(oldArray, newArray, clusterNode->n, oldOffset, newOffset);
      } 
    }
    clusteredPts = newArray;

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

  // tuple<double, bool> getDistNaive(intT i,  intT j){
  //   // double result = FINDNN::dualtree_serial2<nodeT, FAve, double>(getNode(i), getNode(j), fAve);
  //   // return make_tuple(result/getNode(i)->n/getNode(j)->n, false);
  //   return getDistNaive(getNode(i),  getNode(j));
  // }

  ~distAverage1(){
    free(clusteredPts1);
    free(clusteredPts2);
    free(clusterOffsets);
  }

};

// group clustered points at the end of each round
// compute distance by traverse clustered points
template<int dim, class pointTT, class nodeTT>
struct distAverage3: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;

  // private:
  pointT **clusteredPts1;
  pointT **clusteredPts2;//clusteredPts point to one of the two
  pointT *PP;//used to initNodes only
  intT *clusterOffsets;  //same order as activeClusters in finder
  pointT **clusteredPts;

  // copy points from oldArray[oldOffset:oldOffset+n] to newArray[newOffset:newOffset+n]
  inline void copyPoints(pointT **oldArray, pointT **newArray, intT copyn, intT oldOffset, intT newOffset){
      parallel_for(intT j = 0; j < copyn; ++j){
          newArray[newOffset + j] = oldArray[oldOffset+j];
      }
    }

  inline static void printName(){
    cout << "distAverage3 point* array, norm = 1" << endl;
  }
    
  distAverage3(pointT *t_PP, intT n, bool t_no_cache):PP(t_PP){
    clusteredPts1 = newA(pointT *, 2*n);
    clusteredPts2 = clusteredPts1 + n;
    clusteredPts = clusteredPts1;
    parallel_for(intT i=0; i<n; ++i) {
      clusteredPts1[i]=PP + i;
    }
    clusterOffsets = newA(intT, n);
    parallel_for(intT i=0; i<n; ++i) {clusterOffsets[i] = i;}
  }

    inline void initNodes(nodeT *nodes, intT n){
      parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
    }

    double getDistNaive(nodeT *inode,  nodeT *jnode, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    pair<pair<intT, intT>, double> result = FINDNN::bruteForceAverage(inode, jnode, clusteredPts);
    return result.second;
  }

    template<class kdnodeT2, class Fs>
    inline void getRadius(intT cid, kdnodeT2 *root, Fs *fs){
      this->getRadius(cid, root, fs);
    }

  template<class F>
  inline void update(intT round, F *finder){
    intT *activeClusters = finder->activeClusters;
    intT  C = finder->C;
    //  put points into clusters
    parallel_for(intT i = 0; i < C; ++i){
      clusterOffsets[i] = finder->getNode(activeClusters[i])->n;
    }
    sequence::prefixSum(clusterOffsets, 0, C);
    pointT **oldArray = clusteredPts;
    pointT **newArray = clusteredPts1;
    if(clusteredPts == clusteredPts1){
      newArray  = clusteredPts2;
    }

    //TODO: optimize for 2 clusters
    parallel_for(intT i = 0; i < C; ++i){
      intT cid  = activeClusters[i];
      nodeT *clusterNode = finder->getNode(cid);
      intT oldOffset = clusterNode->getOffset(); // must save this before setting new
      intT newOffset = clusterOffsets[i];
      clusterNode->setOffset(newOffset);
      if(clusterNode->round == round){//merged this round, copy left and right
        nodeT *clusterNodeL = clusterNode->left;
        oldOffset = clusterNodeL->getOffset();
        copyPoints(oldArray, newArray, clusterNodeL->n, oldOffset, newOffset);
        newOffset += clusterNodeL->n;

        nodeT *clusterNodeR = clusterNode->right;
        oldOffset = clusterNodeR->getOffset();
        copyPoints(oldArray, newArray, clusterNodeR->n, oldOffset, newOffset);

      }else{//not merged this round, just copy
        copyPoints(oldArray, newArray, clusterNode->n, oldOffset, newOffset);
      } 
    }
    clusteredPts = newArray;

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

  ~distAverage3(){
    free(clusteredPts1);
    // free(clusteredPts2);
    free(clusterOffsets);
  }

};

// compute distance using linkedlist
template<int dim, class pointTT, class nodeTT>
struct distAverage2: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;
  typedef FINDNN::NNaverage<dim, nodeT, pointT *> FAve;

  // private:
  pointT *PP;//used to initNodes only
  FAve *fAve;

  // public:
  LDS::node<pointT *> *PPLList;

  inline static void printName(){
    cout << "distAverage2 point linkedlist, norm = 1" << endl;
  }

  distAverage2(pointT *t_PP, intT n, bool t_no_cache):PP(t_PP){
    PPLList = newA(LDS::node<pointT *>, n);
    parallel_for(intT i=0; i<n; ++i) {PPLList[i] = LDS::node<pointT *>(PP+i, nullptr);}
    fAve = new FAve();
  }

  inline void initNodes(nodeT *nodes, intT n){
    parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PPLList+i);}
  }

  double getDistNaive(nodeT *inode,  nodeT *jnode, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    double result;
    result = FINDNNP::dualtree2<nodeT, FAve, double>(inode, jnode, fAve, result);
    return result/inode->n/jnode->n;
  }

      template<class kdnodeT2, class Fs>
    inline void getRadius(intT cid, kdnodeT2 *root, Fs *fs){
      this->getRadius(cid, root, fs);
    }

  template<class F>
  inline void update(intT round, F *finder){}

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

  ~distAverage2(){
    free(PPLList);
  }

};

template<int dim, class pointTT, class nodeTT, class nodeInfo>
struct distAverage4: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;
  typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
  typedef kdTree<dim, pointT, nodeInfo> kdtreeT;

  bool id_only = true;
  LDS::Method method = LDS::AVGSQ;
  pointT *PP;

  pointT *centers; //ith location is the center of cluster activeCluster[i], need to be contiguous when building tree
  double *vars; //ith location is the variance of cluster cid, sum(x_i-mu)/n
  intT *centerMap; // ith location is the idx of cluster i in centers 
  intT *sizes; //ith location is the size of cluster i

  inline static void printName(){
    cout << "dist Avg 4, use centers and variances, norm = 2" << endl;
  }
    
  distAverage4(pointT *t_PP, intT t_n, bool t_no_cache): PP(t_PP){
    // if(!t_no_cache){
    //   cout << "no cache needs to be true" << endl;
    //   exit(1);
    // }

    centers = newA(pointT, t_n);
    vars = newA(double, t_n);
    sizes = newA(intT, t_n);
    centerMap = newA(intT, t_n);

    parallel_for(intT i=0; i<t_n; ++i) {centers[i] = PP[i];}
    parallel_for(intT i=0; i<t_n; ++i) {vars[i] = 0;}
    parallel_for(intT i=0; i<t_n; ++i) {sizes[i] = 1;}
    parallel_for(intT i=0; i<t_n; ++i) {centerMap[i] = i;}

  }

  inline void initNodes(nodeT *nodes, intT n){
    parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
  }

  inline double getDistNaive(intT cid1, intT cid2, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    intT idx1 = centerMap[cid1];
    intT idx2 = centerMap[cid2];
    return sqrt(centers[idx1].pointDistSq(centers[idx2]) + vars[cid1] + vars[cid2]);
  }

  inline double getDistNaive(nodeT *inode,  nodeT *jnode, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    return getDistNaive(inode->cId, jnode->cId);
  }

  template<class kdnodeT, class Fs>
  void getRadius(intT cid, kdnodeT *root, Fs *fs){
    this->getRadius(cid, root, fs);
  }

  template<class F>
  inline void update(intT round, F *finder){

    //make points array
    intT  C = finder->C;

    parallel_for(intT i = 0; i < C; ++i){
      intT cid = finder->activeClusters[i];
      // intT oldMap = centerMap[cid];
      nodeT *clusterNode = finder->getNode(cid);
      centers[i] = pointT(clusterNode->center, cid);
      centerMap[cid] = i;
      sizes[cid] = clusterNode->n;
      if(clusterNode->round == round){
        vars[cid] = updateVar(clusterNode);
      }
    }

    // build kdtree
    finder->kdtree->kdTreeRebuild(centers, C);
  }

  inline double updateVar(nodeT *clusterNode){
      intT left_id = clusterNode->left->cId;
      intT left_n = clusterNode->left->size();
      pointT left_mu = clusterNode->left->center;
      intT right_id = clusterNode->right->cId;
      intT right_n = clusterNode->right->size();
      pointT right_mu = clusterNode->right->center;
      // ok to use same vars array, only merge two clusters at once, no chaini of merges
      double new_var = left_n * left_mu.pointDistSq(clusterNode->center) + \
                right_n * right_mu.pointDistSq(clusterNode->center) + \
                left_n * vars[left_id] + right_n * vars[right_id]; 
      return new_var/(left_n+right_n);
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
    parallel_for(intT i=0; i<2*finder->n-1; ++i) {finder->nodes[i].height = finder->nodes[i].height * finder->nodes[i].height;}
    parallel_for(intT i=0; i<finder->n; ++i) {finder->uf->values[i] = finder->uf->values[i]  * finder->uf->values[i];}
  }

  ~distAverage4(){
    free(centers);
    free(vars);
    free(sizes);
    free(centerMap);
  }

};

// group clustered points at the end of each round
// compute distance by traverse clustered points
// same as distAverage1, except rebuild kdtree at the end
template<int dim, class pointTT, class nodeTT>
struct distAverage5: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;

  // private:
  pointT *clusteredPts1;
  pointT *clusteredPts2;//clusteredPts point to one of the two
  pointT *PP;//used to initNodes only
  intT *clusterOffsets;  //same order as activeClusters in finder
  pointT *clusteredPts;

  LDS::Method method = LDS::AVG;

  pointT *centers; //ith location is the center of cluster activeCluster[i], need to be contiguous when building tree
  // intT *centerMap; // ith location is the idx of cluster i in centers 

  // copy points from oldArray[oldOffset:oldOffset+n] to newArray[newOffset:newOffset+n]
  inline void copyPoints(pointT *oldArray, pointT *newArray, intT copyn, intT oldOffset, intT newOffset){
      parallel_for(intT j = 0; j < copyn; ++j){
          newArray[newOffset + j] = oldArray[oldOffset+j];
      }
    }
  
  inline static void printName(){
    cout << "distAverage5 point array, norm = 1, rebuild tree" << endl;
  }

  distAverage5(pointT *t_PP, intT n, bool t_no_cache):PP(t_PP){
    clusteredPts1 = newA(pointT, n);
    clusteredPts2 = newA(pointT, n);
    clusteredPts = clusteredPts1;
    parallel_for(intT i=0; i<n; ++i) {
      clusteredPts1[i]=PP[i];
    }
    clusterOffsets = newA(intT, n);
    parallel_for(intT i=0; i<n; ++i) {clusterOffsets[i] = i;}

    centers = newA(pointT, n);
    // centerMap = newA(intT, t_n);

    parallel_for(intT i=0; i<n; ++i) {centers[i] = PP[i];}
    // parallel_for(intT i=0; i<t_n; ++i) {centerMap[i] = i;}

  }

    inline void initNodes(nodeT *nodes, intT n){
      parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
    }

    double getDistNaive(intT cid1, intT cid2, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
        return this->getDistNaive(cid1,cid2);
    }

    double getDistNaive(nodeT *inode,  nodeT *jnode, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    pair<pair<intT, intT>, double> result = FINDNN::bruteForceAverage(inode, jnode, clusteredPts);
    return result.second;
    }
  
    template<class kdnodeT2, class Fs>
    inline void getRadius(intT cid, kdnodeT2 *root, Fs *fs){
      this->getRadius(cid, root, fs);
    }



  template<class F>
  inline void update(intT round, F *finder){
    intT *activeClusters = finder->activeClusters;
    intT  C = finder->C;
    //  put points into clusters
    parallel_for(intT i = 0; i < C; ++i){
      clusterOffsets[i] = finder->getNode(activeClusters[i])->n;
    }
    sequence::prefixSum(clusterOffsets, 0, C);
    pointT *oldArray = clusteredPts;
    pointT *newArray = clusteredPts1;
    if(clusteredPts == clusteredPts1){
      newArray  = clusteredPts2;
    }

    //TODO: optimize for 2 clusters
    parallel_for(intT i = 0; i < C; ++i){
      intT cid  = activeClusters[i];
      nodeT *clusterNode = finder->getNode(cid);
      intT oldOffset = clusterNode->getOffset(); // must save this before setting new
      intT newOffset = clusterOffsets[i];
      clusterNode->setOffset(newOffset);
      if(clusterNode->round == round){//merged this round, copy left and right
        nodeT *clusterNodeL = clusterNode->left;
        oldOffset = clusterNodeL->getOffset();
        copyPoints(oldArray, newArray, clusterNodeL->n, oldOffset, newOffset);
        newOffset += clusterNodeL->n;

        nodeT *clusterNodeR = clusterNode->right;
        oldOffset = clusterNodeR->getOffset();
        copyPoints(oldArray, newArray, clusterNodeR->n, oldOffset, newOffset);

      }else{//not merged this round, just copy
        copyPoints(oldArray, newArray, clusterNode->n, oldOffset, newOffset);
      } 
    }
    clusteredPts = newArray;

    parallel_for(intT i = 0; i < C; ++i){
      intT cid = finder->activeClusters[i];
      nodeT *clusterNode = finder->getNode(cid);
      centers[i] = pointT(clusterNode->center, cid);
      // centerMap[cid] = i;
    }

    // build kdtree
    finder->kdtree->kdTreeRebuild(centers, C);

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
    free(clusteredPts1);
    free(clusteredPts2);
    free(clusterOffsets);
    free(centers);
  }

};

//////////////////////////////////// Complete Linkage ///////////////////////////////

// same as distComplete1, uses array instead of cache
template<int dim, class pointTT, class nodeTT, class nodeInfo>
struct distComplete3: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;
  typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
  typedef kdTree<dim, pointT, nodeInfo> kdtreeT;
  typedef FINDNN::NNcomplete1<dim, kdnodeT> FComp;

  typedef tuple<intT, intT, long> clusterCacheT;//(round, count, cid), use long to pack to 2^i bytes
  // typedef pair<pair<intT, intT>, pair<intT, intT>> clusterCacheT; 

  // private:
  pointT *PP;//used to initNodes only

  kdtreeT **kdtrees;
  clusterCacheT *clusterTbs; // cluster to count
  uintT n;
  intT PNum;
  bool no_cache;
  const bool id_only = true;
  bool nn_process = true;
  LDS::Method method = LDS::COMP;
  intT round = 1;

  inline static void printName(){
    cout << "distComplete3 point* array, kdtree, array for counting" << endl;
  }

  clusterCacheT *initClusterTb(intT pid, intT C){return clusterTbs+(n*pid);}

  struct cmp{
    bool operator () (clusterCacheT i, clusterCacheT j) { // i is new
        if(get<0>(i) == get<0>(j) && get<2>(i) == get<2>(j)) return true;
      return false;
    }
  };

  // template<int k>
  // inline intT get(clusterCacheT &i){
  //   if (k == 0) return i.first.first;
  //   if (k == 1) return i.first.second; 
  //   if (k == 2) return i.second.first;
  // }

  // template<int k>
  // inline void inc(clusterCacheT* tb, intT Rid, intT a){
  //   if (k == 0) tb[Rid].first.first += a; 
  //   else if (k == 1) tb[Rid].first.second += a; 
  //   else if (k == 2) tb[Rid].second.first += a; 
  // }

  // inline clusterCacheT make_entry(intT round, intT cid, intT a){
  //   return make_pair(make_pair(round, a), make_pair(cid, 0));
  // }

  //(count, if reached table end)
  //assume range query is serial
  inline tuple<intT, bool> incrementTable(clusterCacheT* tb, intT Rid, intT cid, intT a = 1){

    // utils::writeMax(tb+Rid, make_tuple(round, 0, (long)cid), cmp());
    // utils::writeAdd( &get<1>(*(tb+Rid)), a);  //((intT *)tb)+(2*Rid + 1)
    if(get<0>(tb[Rid]) != round || get<2>(tb[Rid]) != cid){
      tb[Rid] =  make_tuple(round, a, (long)cid);//make_entry(round, cid, a);
    }else{
      get<1>(tb[Rid]) += a;
    }
    return make_tuple(get<1>(tb[Rid]), false);
  }
    
  distComplete3(pointT *t_PP, intT t_n, bool t_no_cache):PP(t_PP), n(t_n), no_cache(t_no_cache){
    kdtrees = (kdtreeT **) malloc(n*sizeof(kdtreeT *));
    parallel_for(intT i=0; i<n; ++i) {
      kdtrees[i] = new kdtreeT(PP + i, 1, false);
    }

    // if(no_cache){
      PNum =  getWorkers();
      clusterTbs = newA(clusterCacheT, PNum * n);
      parallel_for(intT i = 0; i < PNum * n; ++i){
        clusterTbs[i] = make_tuple(1, 0, (long)-1);
      }
    // }else{
    //   clusterTbs = nullptr;
    // }

  }

    inline void initNodes(nodeT *nodes, intT n){
      parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
    }

  double getDistNaive(intT cid1, intT cid2, double lb = -1, double ub = numeric_limits<double>::max(), bool par = true){ 
    double result;
    if(kdtrees[cid1]->getN() + kdtrees[cid2]->getN() < 200){
      pair<pair<intT, intT>, double> result = bruteForceFarthest(kdtrees[cid1], kdtrees[cid2]);
      return result.second;
    }

    LDS::EDGE e;
    if(lb == -1){
        lb = kdtrees[cid1]->items[0]->pointDist(kdtrees[cid2]->items[0]);
        e = LDS::EDGE(kdtrees[cid1]->items[0]->idx(),kdtrees[cid2]->items[0]->idx(),lb);
        if(lb > ub) return LARGER_THAN_UB;
    }else{
        e = LDS::EDGE(-1,-1,lb);
    }

    FComp fComp = FComp(e, ub);//FComp(LDS::EDGE(cid1,cid2, result));
    if(par){
        FINDNNP::dualtree<kdnodeT, FComp>(kdtrees[cid1]->root, kdtrees[cid2]->root, &fComp); 
    }else{
        FINDNN::dualtree_serial<kdnodeT, FComp>(kdtrees[cid1]->root, kdtrees[cid2]->root, &fComp); 
    }    
    result = fComp.getResultW();
    return result;
  }

   double getDistNaive(nodeT *inode,  nodeT *jnode, double lb = -1, double ub = numeric_limits<double>::max(), bool par = true){ 
    return getDistNaive(inode->cId, jnode->cId, lb, ub, par);
  }

  template<class kdnodeT, class Fs>
  void getRadius(intT cid, kdnodeT *root, Fs *fs){
    FINDNNP::dualtree<kdnodeT, Fs>(kdtrees[cid]->root, root, fs); 
  }

  template<class F>
  inline void update(intT _round, F *finder){
    round = _round;
    intT *activeClusters = finder->activeClusters;
    intT  C = finder->C;

    parallel_for(intT i = 0; i < C; ++i){
      intT cid  = activeClusters[i];
      nodeT *clusterNode = finder->getNode(cid);
      if(clusterNode->round == round){//merged this round, build new tree
        intT cid1 = clusterNode->left->cId;
        intT cid2 = clusterNode->right->cId;
        intT n1 = clusterNode->left->n;
        intT n2 = clusterNode->right->n;
        intT newN = n1 + n2;
        pointT **items = newA(pointT *, newN);
        parallel_for(intT i =0; i < n1; ++i){
          items[i] = kdtrees[cid1]->items[i];
        }
        parallel_for(intT i =0; i < n2; ++i){
          items[n1+i] = kdtrees[cid2]->items[i];
        }
        
        delete kdtrees[cid1];
        delete kdtrees[cid2];
        if(C > 1) {kdtrees[cid] = new kdTree<dim, pointT, nodeInfo>(items, newN, newN > 2000);}
        else{free(items);}
      }
      
    }
    round++;// when using it, we want to use the next round
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
    free(kdtrees);
    // if(no_cache){
      free(clusterTbs);
    // }
  }

};



//////////////////////////////////// Ward's Linkage ///////////////////////////////

// pointer array
// not keeping clustered points together in 1 array, 
// because when rebuidling kdtrees, need to allocate nodes anyway
// TODO: arrays to copy both nodesm and points
// switch by looking at how many merged clusters
// build kdtrees of center points each round
template<int dim, class pointTT, class nodeTT, class nodeInfo, class M>
struct distWard1: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;
  typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
  typedef kdTree<dim, pointT, nodeInfo> kdtreeT;

  bool no_cache = true;
  bool id_only = true;
  LDS::Method method = LDS::WARD;
  pointT *PP;
  M marker;

  // kdtreeT *kdtree;
  pointT *centers; //ith location is the center of cluster activeCluster[i]
  intT *centerMap; // ith location is the idx of cluster i in centers
  intT *sizes; //ith location is the size of cluster i


  inline static void printName(){
    cout << "dist Ward 1" << endl;
  }
    
  distWard1(pointT *t_PP, intT t_n, bool t_no_cache): PP(t_PP){
    // if(!no_cache){
    //   cout << "no cache needs to be true" << endl;
    //   exit(1);
    // }
    centers = newA(pointT, t_n);
    sizes = newA(intT, t_n);
    centerMap = newA(intT, t_n);

    parallel_for(intT i=0; i<t_n; ++i) {centers[i] = PP[i];}
    parallel_for(intT i=0; i<t_n; ++i) {sizes[i] = 1;}
    parallel_for(intT i=0; i<t_n; ++i) {centerMap[i] = i;}


    // kdtree = new kdtreeT(centers, t_n);
    marker = M(sizes);

  }

  inline void initNodes(nodeT *nodes, intT n){
    parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
  }

  inline double getDistNaive(intT cid1, intT cid2, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    // return this->getDistNaive(cid1,cid2);
    double ni = (double) sizes[cid1]; 
    double nj = (double) sizes[cid2];
    // if(ni + nj > 2) 
    return sqrt(2*(ni*nj)*centers[centerMap[cid1]].pointDistSq(centers[centerMap[cid2]])/(ni + nj));
    // return centers[centerMap[cid1]].pointDist(centers[centerMap[cid2]]);
  }

  inline double getDistNaive(nodeT *inode,  nodeT *jnode, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    double ni = (double) inode->n; 
    double nj = (double) jnode->n;
    if(ni + nj > 2) return sqrt(2*(ni*nj)*inode->dist(jnode)/(ni + nj));
    return sqrt(inode->dist(jnode));
  }

  template<class kdnodeT, class Fs>
  void getRadius(intT cid, kdnodeT *root, Fs *fs){
    this->getRadius(cid, root, fs);
  }

  template<class F>
  inline void update(intT round, F *finder){
    //make points array
    intT  C = finder->C;
    
    parallel_for(intT i = 0; i < C; ++i){
      intT cid = finder->activeClusters[i];
      centers[i] = pointT(finder->getNode(cid)->center, cid);
      centerMap[cid] = i;
      sizes[cid] = finder->getNode(cid)->n;
    }

    // build kdtree
    finder->kdtree->kdTreeRebuild(centers, C);

    //mark min sizes
    FINDNN::singletree<kdnodeT, M, typename M::infoT>(finder->kdtree->root, &marker, marker.initVal);
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
    // delete kdtree;
    free(centers);
    free(sizes);
    free(centerMap);
  }

};


// pointer array, dynamic kdtree of centers
template<int dim, class pointTT, class nodeTT, class nodeInfo, class M>
struct distWard2: public distAbstract {
  typedef pointTT pointT;
  typedef nodeTT nodeT;
  typedef dynamicKdNode<dim, pointT, nodeInfo> kdnodeT;
  typedef dynamicKdTree<dim, pointT, nodeInfo> kdtreeT;

  bool no_cache = true;
  bool id_only = true;
  pointT *PP;
  M marker;

  pointT *centers; //ith location is the center of cluster activeCluster[i]
  pointT *centersDepr; //(2i,2i+1)th location is the center of cluster activeCluster[i]
  bool *flag;
  intT *centerMap; // ith location is the idx of cluster i in centers
  intT *sizes; //ith location is the size of cluster i
  pointT *treePts;
  intT counter = 0;


  inline static void printName(){
    cout << "dist Ward 2 dynamic kdtree" << endl;
  }
    
  distWard2(pointT *t_PP, intT t_n, bool t_no_cache): PP(t_PP){
    // if(!no_cache){
    //   cout << "no cache needs to be true" << endl;
    //   exit(1);
    // }
    centers = newA(pointT, t_n);
    centersDepr = newA(pointT, 2*t_n);
    treePts = newA(pointT, 2*t_n);
    flag = newA(bool, 3*t_n);
    sizes = newA(intT, t_n);
    centerMap = newA(intT, t_n);

    parallel_for(intT i=0; i<t_n; ++i) {centers[i] = PP[i];}
    parallel_for(intT i=0; i<t_n; ++i) {sizes[i] = 1;}
    parallel_for(intT i=0; i<t_n; ++i) {centerMap[i] = i;}


    // kdtree = new kdtreeT(centers, t_n);
    marker = M(sizes);

  }

  inline void initNodes(nodeT *nodes, intT n){
    parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
  }

  inline double getDistNaive(intT cid1, intT cid2, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    // return this->getDistNaive(cid1,cid2);
    double ni = (double) sizes[cid1]; 
    double nj = (double) sizes[cid2];
    // if(ni + nj > 2) 
    return sqrt(2*(ni*nj)*centers[centerMap[cid1]].pointDistSq(centers[centerMap[cid2]])/(ni + nj));
    // return centers[centerMap[cid1]].pointDist(centers[centerMap[cid2]]);
  }

  inline double getDistNaive(nodeT *inode,  nodeT *jnode, 
                          double lb = -1, double ub = numeric_limits<double>::max(), 
                          bool par = true){ 
    double ni = (double) inode->n; 
    double nj = (double) jnode->n;
    if(ni + nj > 2) return sqrt(2*(ni*nj)*inode->dist(jnode)/(ni + nj));
    return sqrt(inode->dist(jnode));
  }

  template<class kdnodeT, class Fs>
  void getRadius(intT cid, kdnodeT *root, Fs *fs){
    this->getRadius(cid, root, fs);
  }

  template<class F>
  inline void update(intT round, F *finder){
    //make points array
    intT  C = finder->C;

    parallel_for(intT i = 0; i < C; ++i){
      intT cid  = finder->activeClusters[i];
      nodeT *clusterNode = finder->getNode(cid);
      centerMap[cid] = i;
      centers[i] = pointT(clusterNode->center, cid);
      sizes[cid] = finder->getNode(cid)->n;
      if(clusterNode->round == round){//merged this round, build new tree
        centersDepr[2*i] = pointT(clusterNode->left->center, clusterNode->left->cId);
        centersDepr[2*i+1] = pointT(clusterNode->right->center, clusterNode->right->cId);
        flag[i] = true;
        flag[C+2*i] = true;
        flag[C+2*i+1] = true;
      }else{
        centersDepr[2*i] = pointT(clusterNode->center, -1);
        centersDepr[2*i+1] = pointT(clusterNode->center, -1);
        flag[i] = false;
        flag[C+2*i] = false;
        flag[C+2*i+1] = false;
      }
    }    
    _seq<pointT> addPts = sequence::pack<pointT, intT, sequence::getA<pointT,intT> >(treePts+counter, flag, 0, C, sequence::getA<pointT,intT>(centers));
    _seq<pointT> removePts = sequence::pack<pointT, intT, sequence::getA<pointT,intT> >(NULL, flag+C, 0, 2*C, sequence::getA<pointT,intT>(centersDepr));

    counter += addPts.n;
    if(addPts.n > C/2){
      // build kdtree
      delete finder->kdtree;
      finder->kdtree = new kdtreeT(centers, C);//Don't pass in PP! we need centers, not points
      counter = 0;
    }else{
      finder->kdtree->erase(removePts.A, removePts.n);
      finder->kdtree->insert(addPts.A, addPts.n);
    }


    // free(addPts.A);
    free(removePts.A);

    //mark min sizes
    FINDNN::singletree<kdnodeT, M, typename M::infoT>(finder->kdtree->root, &marker, marker.initVal);
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



  ~distWard2(){
    // delete kdtree;
    free(centers);
    free(centersDepr);
    free(sizes);
    free(flag);
    free(centerMap);
    free(treePts);
  }

};

//////////////////////// HDBSCAN - Boruvka //////////////////////////////////

// default k to 2
template<int dim, class pointTT>
struct distHDBSCAN1{
  typedef pointTT pointT;

  // private:
  pointT *PP;//used to initNodes only
  double *coredists;
  int k;

  inline static void printName(){
    cout << "distHDBSCAN1 point array" << endl;
  }

  distHDBSCAN1(pointT *t_PP, intT n, commandLine params):PP(t_PP){
    k = params.getOptionIntValue("-k",2);
  }

  template<class kdtreeT>
  inline void initTree(kdtreeT* kdtree, intT n){
    pointT** A = newA(pointT*, k*n);
    coredists = newA(double, n);
    parallel_for(intT i = 0; i < n; ++i){
        kdtree->kNN(&PP[i], k, A+i*k);
        coredists[i] = PP[i].pointDist(A[(i+1)*k-1]);
    }
    free(A);
    using M = FINDNN::MarkKdTreeCoreDist<dim, typename kdtreeT::nodeT>;
    M marker = M(coredists);
    FINDNN::singletree<typename kdtreeT::nodeT, M, typename M::infoT>(kdtree->root, &marker, marker.initVal);
  }

  inline double getPointDist(const intT &i, const intT &j){
    double dist = PP[i].pointDist(PP[j]);
    return max(dist, max(coredists[i], coredists[j]));
  }

  inline double getPointDist(const pointT &p, const pointT &q){
    double dist = p.pointDist(q);
    return max(dist, max(coredists[p.i], coredists[q.i]));
  }

  inline double getPointDist(pointT *p, pointT *q){
    double dist = p->pointDist(q);
    return max(dist, max(coredists[p->i], coredists[q->i]));
  }

  inline bool skipBaseCase(const double &ub, pointT *q){
    return coredists[q->i] > ub; // q can't update the nn of cluster(q)
  }
  
  inline bool useSimpleNN(intT round){
    return false;
  }

  template<class kdnodeT>
  inline void updateUB(kdnodeT *Q){
    // if(Q->nInfo.max_core_dist == -1){
    //   if(Q->isLeaf()){
    //     double temp = coredists[Q->items[0]->idx()];
    //     for(intT i=1; i<Q->size(); ++i){
    //       double q_cd = coredists[Q->items[i]->idx()]; 
		//       if(q_cd > temp) temp = q_cd;
    //     }
    //      Q->nInfo.max_core_dist = temp;
    //   }else{
    //     Q->nInfo.max_core_dist = max((Q->left->nInfo).max_core_dist, (Q->right->nInfo).max_core_dist);
    //   }
    // }//end of max_core_dist == -1
    if ( Q->nInfo.getCId() < 0){
        (Q->nInfo).updateUB(max(Q->Diameter(),Q->nInfo.max_core_dist) );
    }
  }

  ~distHDBSCAN1(){
    free(coredists);
  }

};

template<int dim, class pointTT>
struct distSingle {
  typedef pointTT pointT;

  // private:
  pointT *PP;//used to initNodes only

  inline static void printName(){
    cout << "distSingle" << endl;
  }

  distSingle(pointT *t_PP, intT n, commandLine params):PP(t_PP){
  }

  template<class kdtreeT>
  inline void initTree(kdtreeT* kdtree, intT n){
  }

  inline double getPointDist(const intT &i, const intT &j){
    return PP[i].pointDist(PP[j]);
  }

  inline double getPointDist(const pointT &p, const pointT &q){
    return p.pointDist(q);
  }

  inline double getPointDist(pointT *p, pointT *q){
    return p->pointDist(q);
  }

  inline bool skipBaseCase(double ub, pointT *q){
    return false;
  }

  inline bool useSimpleNN(intT round){
    return round == 1;
  }

  template<class kdnodeT>
  inline void updateUB(kdnodeT *Q){
    if ( Q->nInfo.getCId() < 0){
      (Q->nInfo).updateUB(Q->Diameter());
    }
  }

  ~distSingle(){ }

};


#endif