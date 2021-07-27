#ifdef PERF_RANGE
  double range_time = 0;
  double find_nns_time = 0;
  // intT cand_num = 0;
  double max_round_time = 0;
  double max_range_time = 0;
  double max_nns_time = 0;

inline void report_perf_range(intT round){
    if(round > 1 && (round < 5 || round % PRINT_FREQ == 0)){
      UTIL::PrintFunctionItem("PERF_RANGE", "nns", find_nns_time);
      UTIL::PrintFunctionItem("PERF_RANGE", "init_tree_time", init_tree_time);
      // UTIL::PrintFunctionItem("PERF_RANGE", "tree_bb_time", tree_bb_time);
      UTIL::PrintFunctionItem("PERF_RANGE", "range+nnc", range_time);
      UTIL::PrintFunctionItem("PERF_RANGE", "cand_num", cand_num);
  }
}

inline void reset_perf_range(){
    range_time = 0;
    find_nns_time = 0;
    init_tree_time =  0;
    // tree_bb_time = 0;
    cand_num = 0;
  }
#endif


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


// // pointer array, dynamic kdtree of centers
// template<int dim, class pointTT, class nodeTT, class nodeInfo, class M>
// struct distWard2: public distAbstract {
//   typedef pointTT pointT;
//   typedef nodeTT nodeT;
//   typedef dynamicKdNode<dim, pointT, nodeInfo> kdnodeT;
//   typedef dynamicKdTree<dim, pointT, nodeInfo> kdtreeT;

//   bool no_cache = true;
//   bool id_only = true;
//   pointT *PP;
//   M marker;

//   pointT *centers; //ith location is the center of cluster activeCluster[i]
//   pointT *centersDepr; //(2i,2i+1)th location is the center of cluster activeCluster[i]
//   bool *flag;
//   intT *centerMap; // ith location is the idx of cluster i in centers
//   intT *sizes; //ith location is the size of cluster i
//   pointT *treePts;
//   intT counter = 0;


//   inline static void printName(){
//     cout << "dist Ward 2 dynamic kdtree" << endl;
//   }
    
//   distWard2(pointT *t_PP, intT t_n, bool t_no_cache): PP(t_PP){
//     // if(!no_cache){
//     //   cout << "no cache needs to be true" << endl;
//     //   exit(1);
//     // }
//     centers = newA(pointT, t_n);
//     centersDepr = newA(pointT, 2*t_n);
//     treePts = newA(pointT, 2*t_n);
//     flag = newA(bool, 3*t_n);
//     sizes = newA(intT, t_n);
//     centerMap = newA(intT, t_n);

//     parallel_for(intT i=0; i<t_n; ++i) {centers[i] = PP[i];}
//     parallel_for(intT i=0; i<t_n; ++i) {sizes[i] = 1;}
//     parallel_for(intT i=0; i<t_n; ++i) {centerMap[i] = i;}


//     // kdtree = new kdtreeT(centers, t_n);
//     marker = M(sizes);

//   }

//   inline void initNodes(nodeT *nodes, intT n){
//     parallel_for(intT i=0; i<n; ++i) {nodes[i] = nodeT(i, PP[i]);}
//   }

//   inline double getDistNaive(intT cid1, intT cid2, 
//                           double lb = -1, double ub = numeric_limits<double>::max(), 
//                           bool par = true){ 
//     // return this->getDistNaive(cid1,cid2);
//     double ni = (double) sizes[cid1]; 
//     double nj = (double) sizes[cid2];
//     // if(ni + nj > 2) 
//     return sqrt(2*(ni*nj)*centers[centerMap[cid1]].pointDistSq(centers[centerMap[cid2]])/(ni + nj));
//     // return centers[centerMap[cid1]].pointDist(centers[centerMap[cid2]]);
//   }

//   inline double getDistNaive(nodeT *inode,  nodeT *jnode, 
//                           double lb = -1, double ub = numeric_limits<double>::max(), 
//                           bool par = true){ 
//     double ni = (double) inode->n; 
//     double nj = (double) jnode->n;
//     if(ni + nj > 2) return sqrt(2*(ni*nj)*inode->dist(jnode)/(ni + nj));
//     return sqrt(inode->dist(jnode));
//   }

//   template<class kdnodeT, class Fs>
//   void getRadius(intT cid, kdnodeT *root, Fs *fs){
//     this->getRadius(cid, root, fs);
//   }

//   template<class F>
//   inline void update(intT round, F *finder){
//     //make points array
//     intT  C = finder->C;

//     parallel_for(intT i = 0; i < C; ++i){
//       intT cid  = finder->activeClusters[i];
//       nodeT *clusterNode = finder->getNode(cid);
//       centerMap[cid] = i;
//       centers[i] = pointT(clusterNode->center, cid);
//       sizes[cid] = finder->getNode(cid)->n;
//       if(clusterNode->round == round){//merged this round, build new tree
//         centersDepr[2*i] = pointT(clusterNode->left->center, clusterNode->left->cId);
//         centersDepr[2*i+1] = pointT(clusterNode->right->center, clusterNode->right->cId);
//         flag[i] = true;
//         flag[C+2*i] = true;
//         flag[C+2*i+1] = true;
//       }else{
//         centersDepr[2*i] = pointT(clusterNode->center, -1);
//         centersDepr[2*i+1] = pointT(clusterNode->center, -1);
//         flag[i] = false;
//         flag[C+2*i] = false;
//         flag[C+2*i+1] = false;
//       }
//     }    
//     _seq<pointT> addPts = sequence::pack<pointT, intT, sequence::getA<pointT,intT> >(treePts+counter, flag, 0, C, sequence::getA<pointT,intT>(centers));
//     _seq<pointT> removePts = sequence::pack<pointT, intT, sequence::getA<pointT,intT> >(NULL, flag+C, 0, 2*C, sequence::getA<pointT,intT>(centersDepr));

//     counter += addPts.n;
//     if(addPts.n > C/2){
//       // build kdtree
//       delete finder->kdtree;
//       finder->kdtree = new kdtreeT(centers, C);//Don't pass in PP! we need centers, not points
//       counter = 0;
//     }else{
//       finder->kdtree->erase(removePts.A, removePts.n);
//       finder->kdtree->insert(addPts.A, addPts.n);
//     }


//     // free(addPts.A);
//     free(removePts.A);

//     //mark min sizes
//     FINDNN::singletree<kdnodeT, M, typename M::infoT>(finder->kdtree->root, &marker, marker.initVal);
//   }

//   inline double updateDistO(double dik, double djk, double ni, double nj, double nk, double dij){
//     double ntotal = ni + nj + nk;
//     double d = sqrt( ( ((ni + nk)  * dik * dik) + ((nj + nk) * djk * djk) - (nk * dij * dij) )/ ntotal );
//     return d;
//   }

//   inline double updateDistN(double dikl, double dikr, double djkl, double djkr, 
//                        double ni, double nj, double nkl, double nkr,
//                        double dij, double dklr ){
//     double dik = updateDistO(dikl, dikr, nkl, nkr, ni, dklr);
//     double djk = updateDistO(djkl, djkr, nkl, nkr, nj, dklr);
//     return updateDistO(dik, djk, ni, nj, (nkl + nkr), dij);
//   }



//   ~distWard2(){
//     // delete kdtree;
//     free(centers);
//     free(centersDepr);
//     free(sizes);
//     free(flag);
//     free(centerMap);
//     free(treePts);
//   }

// };
