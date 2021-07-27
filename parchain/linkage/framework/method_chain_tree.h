#ifndef METHOD_CHAIN_TREE_AVE_H
#define METHOD_CHAIN_TREE_AVE_H

#include "sequence.h"
#include "geometry.h"
#include "sampleSort.h"
#include "gettime.h"
#include "unionfind.h"

// #define A_HASH_LINKAGE_PROB_PRINT

#include "ndHash.h"
#include "serialHash.h"

#include "kdTree2.h"

// #define TIMING2
#define ALLOCDOUBLE
#define MAX_CACHE_TABLE_SIZE_INIT 64
#define LINKAGE_LOADFACTOR 2.0
// #define INTT_MAX numeric_limits<intT>::max()
// #define CHECK_NO_CACHE(i) if(no_cache){cout << "no cache " << i << endl; exit(1);}
#define CHECK_NO_CACHE(i)

#include <limits>

#include "shared.h"
#include "chain.h"
#include "cluster.h"
#include "dist.h"
#include "neighbor.h"
#include "neighbor_parallel.h"
#include "rangeQuery.h"

// #define PRINT_FREQ 1
// #define LINKAGE_DOPRINT(round) (round < 5 || round % PRINT_FREQ == 0)

using namespace std;

//TODO: add nodoeIdx to Info 
template<int dim, class distF, class Fr, class M>
class TreeNNFinder { //: public NNFinder<dim>
  public:
  // typedef iPoint<dim> pointT;
  typedef typename distF::pointT pointT;
  // typedef typename FINDNN::CLinkNodeInfo nodeInfo;
  typedef typename Fr::nodeInfo nodeInfo;
  // typedef FINDNN::AveLinkNodeInfo<dim, pointT *> nodeT;
  typedef typename distF::nodeT nodeT;
  // typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
  // typedef kdTree<dim, pointT, nodeInfo > treeT;
  typedef typename Fr::kdnodeT kdnodeT;
  typedef typename Fr::kdtreeT treeT;
  typedef LDS::distCacheT distCacheT;
  typedef LDS::hashClusterAveET hashClusterAveET;

  bool no_cache;
  intT C;// number of clusters
  intT n;// number of points
  pointT *PP;
  intT *rootIdx;
  UnionFind::ParUF<intT> *uf;
  intT *activeClusters; // ids of connected components
  bool *flag;//used in updateActivateClusters
  treeT *kdtree;
  LDS::EDGE *edges; //edges[i] stores the min neigh of cluster i

  distCacheT **cacheTbs; // distance to clusters, store two copies of each distance
#ifdef ALLOCDOUBLE
  uintT hashTableSize=0;
  distCacheT::eType *TA;
#endif

  nodeT *nodes; // preallocated space to store tree nodes
  LDS::edgeComparator2 EC2;
  double eps;
  atomic<intT> nodeIdx; // the index of next node to use for cluster trees

  distF *distComputer;
  M marker;

  intT NAIVE_THRESHOLD = 50;
  intT MAX_CACHE_TABLE_SIZE = 128;
  // bool fullTableFlag = false;


#ifdef BENCHCACHE
  // intT alloc_cache=0;
  intT used_cache=0;
  // intT *allocCacheCounters;
  intT *usedCacheCounters;
#endif
#ifdef PERF_RANGE
  long *distance_computed;
  long *pointsInRange;
  // long max_max_cache_used_allrounds=0;
  // long max_avg_cache_used_allrounds = 0;
  // long times_cache_full;
#endif  

  //TODO: batch allocate cache table space
  //TODO: free cache table when merges
  TreeNNFinder(intT t_n, point<dim>* t_P, UnionFind::ParUF<intT> *t_uf, 
    bool t_noCache, double t_eps, intT t_naive_thresh, intT t_cache_size): 
    uf(t_uf), n(t_n), no_cache(t_noCache), eps(t_eps), NAIVE_THRESHOLD(t_naive_thresh), MAX_CACHE_TABLE_SIZE(t_cache_size){
    EC2 = LDS::edgeComparator2(eps);
    C = n;
    activeClusters = newA(intT, n);
    parallel_for(intT i=0; i<n; ++i) {activeClusters[i] = i;}
    PP = newA(pointT, n);
    parallel_for(intT i=0; i<n; ++i) {PP[i] = pointT(t_P[i], i);}

    distComputer = new distF(PP, n, no_cache);

    rootIdx = newA(intT, n);
    parallel_for(intT i=0; i<n; ++i) {rootIdx[i] = i;}
    nodes = newA(nodeT, 2*n);
    distComputer->initNodes(nodes, n);
    
    flag = newA(bool, C);

    edges = (LDS::EDGE *)malloc(sizeof(LDS::EDGE) * n);
    parallel_for(intT i=0; i<n; ++i) {edges[i] = LDS::EDGE(-1,-1,numeric_limits<double>::max());}

#ifdef ALLOCDOUBLE
    if(!no_cache){
    hashTableSize = min(n, 1 << utils::log2Up((uintT)(LINKAGE_LOADFACTOR*(uintT)MAX_CACHE_TABLE_SIZE)));
    TA = newA(distCacheT::eType, (long)n*2*hashTableSize);
    distCacheT::eType emptyval = LDS::hashClusterAve().empty();
    cacheTbs = newA(distCacheT *, 2*n);
    parallel_for(intT i=0; i<n; ++i) {
      cacheTbs[i] = new distCacheT(min(MAX_CACHE_TABLE_SIZE, min(MAX_CACHE_TABLE_SIZE_INIT, C)), TA + (long)i*hashTableSize, LDS::hashClusterAve(), LINKAGE_LOADFACTOR, true);
    }
    parallel_for(intT i=n; i<2*n; ++i) {
      cacheTbs[i] = new distCacheT(min(MAX_CACHE_TABLE_SIZE, C), TA + (long)i*hashTableSize, LDS::hashClusterAve(), LINKAGE_LOADFACTOR);
    }
    parallel_for(long i=(long)n*hashTableSize; i<(long)n*2*hashTableSize; ++i) {
      TA[i] = emptyval;
    }
#else
    if(!no_cache){
    cacheTbs = newA(distCacheT *, 2*n);
    parallel_for(intT i=0; i<n; ++i) {
      cacheTbs[i] = new distCacheT(min(MAX_CACHE_TABLE_SIZE_INIT, C), LDS::hashClusterAve(), LINKAGE_LOADFACTOR);
    }
#endif

#ifdef BENCHCACHE
    // alloc_cache=n*cacheTbs[0]->m;
    used_cache=0; 
    // allocCacheCounters=newA(intT,n);
    // parallel_for(intT i=0; i<n; ++i) {allocCacheCounters[i]=min(MAX_CACHE_TABLE_SIZE_INIT, C);}
    usedCacheCounters=newA(intT,n);
    parallel_for(intT i=0; i<n; ++i) {usedCacheCounters[i]=0;}
    cout << "alloc cache: " << n*2*hashTableSize << endl;
    cout << "used cache: " << 0 << endl;
#endif
    } // end of if ! no cache

#ifdef PERF_RANGE
  cout << ELTPERCACHELINE * getWorkers() << endl;
  distance_computed = newA(long, ELTPERCACHELINE * getWorkers());
  pointsInRange = newA(long, ELTPERCACHELINE * getWorkers());
  for(intT i=0; i<ELTPERCACHELINE * getWorkers(); ++i) {distance_computed[i]=0;pointsInRange[i]=0;}
#endif  

    nodeIdx.store(n); // have used n nodes
    kdtree = new treeT(PP, n);
    marker = M(uf, nodes, rootIdx);

  }

  ~TreeNNFinder(){
    if(!no_cache){
#ifndef ALLOCDOUBLE
      cacheTbs[rootIdx[activeClusters[0]]]->del();
      delete cacheTbs[rootIdx[activeClusters[0]]];
      cacheTbs[rootIdx[activeClusters[0]]] = nullptr;
      // parallel_for(intT i=0; i<2*n; ++i) {
      //   if(cacheTbs[i] != nullptr) cout << "cache table not freed " << i << endl;
      // }
#else
      parallel_for(intT i=0; i<2*n; ++i) {
        delete cacheTbs[i];
      }
      free(TA);
#endif
      free(cacheTbs);
#ifdef BENCHCACHE
    free(usedCacheCounters);
#endif
    }
    free(PP);
    free(flag);
    free(activeClusters);
    delete kdtree;
    free(nodes);
    free(edges);
    free(rootIdx);
    delete distComputer;
#ifdef PERF_RANGE
    free(distance_computed);
    free(pointsInRange);
#endif 
  }

  inline intT cid(nodeT* node){ return node->cId;}
  inline intT idx(nodeT* node){ return node->idx;}
  inline nodeT *getNode(intT cid){return nodes+rootIdx[cid];}
  inline distCacheT *getTable(intT cid){return cacheTbs[rootIdx[cid]];}
  inline intT idx(intT cid){return idx(getNode(cid));}
  inline intT leftIdx(intT cid){return getNode(cid)->left->idx;}
  inline intT rightIdx(intT cid){return getNode(cid)->right->idx;}
  inline intT cid(intT idx){return cid(nodes[idx]);}


  inline double find(intT qid, intT rid, intT qIdx = -1, intT  rIdx = -1){
    CHECK_NO_CACHE(259)
    if(qIdx == -1){
      qIdx = idx(qid);
      rIdx = idx(rid);
    }

    typename distCacheT::eType result;
    bool reach_thresh;
    tie(result, reach_thresh) = cacheTbs[qIdx]->find_thresh(rid);
    if(!reach_thresh && result.idx == rIdx){
	    return result.dist;
    }
    
    tie(result, reach_thresh) = cacheTbs[rIdx]->find_thresh(qid);
    if(!reach_thresh && result.idx == qIdx){
	    return result.dist;
    }
    return UNFOUND_TOKEN;
  }

  inline double find(nodeT *inode, nodeT *jnode){
	  return find(cid(inode), cid(jnode), idx(inode), idx(jnode));
  }


  inline void insert_helper(intT qid, intT rid, double d, distCacheT *tb1, distCacheT *tb2){
    if(d == LARGER_THAN_UB){
        return;
    }
    CHECK_NO_CACHE(284)
    intT qIdx = idx(qid);
    intT rIdx = idx(rid);
    tb1->insert2(hashClusterAveET(rid, rIdx, d));
    tb2->insert2(hashClusterAveET(qid, qIdx, d));
  }

  inline void insert(intT qid, intT rid, double d){
    if(d == LARGER_THAN_UB){
        return;
    }
    CHECK_NO_CACHE(291)
    intT qIdx = idx(qid);
    intT rIdx = idx(rid);

    cacheTbs[qIdx]->insert2(hashClusterAveET(rid, rIdx, d));
    cacheTbs[rIdx]->insert2(hashClusterAveET(qid, qIdx, d));
  }

  inline double getDistNaive(nodeT *inode,  nodeT *jnode, double lb = -1, double ub = numeric_limits<double>::max(), bool par = true){
    return distComputer->getDistNaive(inode, jnode, lb, ub, par);
  }

  inline double getDistNaive(intT i, intT j, double lb = -1, double ub = numeric_limits<double>::max(), bool par = true){
    return distComputer->getDistNaive(i, j, lb, ub, par);
  }

  inline double getDistNaiveDebug(intT i, intT j){
    return distComputer->getDistNaive(getNode(i), getNode(j), -1,numeric_limits<double>::max(), false );
  }

  // i, j are cluster ids
  // find distance in table 
  // if not found, compute by bruteforce and insert 
  // true if found in table
  // used in merging 
  // range search has its own updateDist
  tuple<double, bool> getDist(intT i,  intT j, double lb = -1, double ub = numeric_limits<double>::max(), bool par = true){
    if(!no_cache){
    double d = find(i, j);
    // if(d == CHECK_TOKEN){cout << "find check token" << endl; exit(1);} // might find is in singleNN step in getNN
    if(d != UNFOUND_TOKEN && d != CHECK_TOKEN) return make_tuple(d, true);
    }
    if(distComputer->id_only) return make_tuple(getDistNaive(i,  j, lb, ub, par), false);
    return make_tuple(getDistNaive(getNode(i),  getNode(j), lb, ub, par), false);
  }

  tuple<double, bool> getDist(nodeT *inode,  nodeT *jnode, double lb = -1, double ub = numeric_limits<double>::max(), bool par = true){
    if(!no_cache){
    double d = find(inode, jnode);
    // if(d == CHECK_TOKEN){cout << "find check token" << endl; exit(1);} // might find is in singleNN step in getNN
    if(d != UNFOUND_TOKEN && d != CHECK_TOKEN) return make_tuple(d, true); 
    }
    return make_tuple(getDistNaive(inode, jnode, lb, ub, par), false);
  }

  //newc is a newly merged cluster
  // newc is new, rid is old
  inline double getNewDistO(intT newc, intT rid){
    nodeT* ql = getNode(newc)->left;
    intT nql = ql->n;
    nodeT* qr = getNode(newc)->right;
    intT nqr = qr->n;
    nodeT* rroot = getNode(rid);
    intT nr = rroot->n;
    double dij = getNode(newc)->getHeight();

    // double n1 = (double)nql * (double)nr;
    // double n2 = (double)nqr * (double)nr;
    // double alln = n1 + n2 ;

    double d1,d2; bool intable;
    tie(d1, intable) = getDist(ql,rroot);
    // d1 = n1 / alln * d1;

    tie(d2, intable) = getDist(qr,rroot);
    // d2 = n2 / alln * d2;

    // return d1 + d2;
    return distComputer->updateDistO(d1, d2, nql, nqr, nr, dij);
  }


  //newc is a newly merged cluster
  // newc is new, rid is merged
  //todo: cids on cluster tree do not need to be marked
  // overwrite entry?
  // TODO: consider changing hashtable to idx -> (idx, dist)
  inline double getNewDistN(intT newc, intT rid){
    nodeT* ql = getNode(newc)->left;
    intT nql = ql->n;
    nodeT* qr = getNode(newc)->right;
    intT nqr = qr->n;
    double dij = getNode(newc)->getHeight();

    nodeT* rl = getNode(rid)->left;
    intT nrl = rl->n;
    nodeT* rr = getNode(rid)->right;
    intT nrr = rr->n;
    double dklr = getNode(rid)->getHeight();

    double d1,d2, d3, d4; bool intable;
    tie(d1, intable) = getDist(ql,rl);
    // d1 = n1 / alln * d1;

    tie(d2, intable) = getDist(ql,rr);
    // d2 = n2 / alln * d2;

    tie(d3, intable) = getDist(qr,rl);
    // d3 = n3 / alln * d3;

    tie(d4, intable) = getDist(qr,rr);
    // d4 = n4 / alln * d4;

    // return d1  + d2  + d3 + d4;
    return distComputer->updateDistN(d1, d2, d3, d4, nql, nqr, nrl, nrr, dij, dklr);
  }

  // store the closest nn in edges
  // assume edges[cid] already has max written
  void getNN_naive(intT cid, double ub = numeric_limits<double>::max(), intT t_nn = -1){
    utils::writeMin(&edges[cid], LDS::EDGE(cid, t_nn, ub), EC2); 
    parallel_for(intT i = 0; i < C; ++i){
      intT cid2 = activeClusters[i];
      if(cid2 != cid){//if(cid2 < cid) { // the larger one might not be a terminal node only work if C == |terminal node|
        double tmpD;
        bool intable;
        tie(tmpD, intable) = getDist(cid, cid2, -1, edges[cid].getW(), true);
        utils::writeMin(&edges[cid], LDS::EDGE(cid,cid2,tmpD), EC2);
        utils::writeMin(&edges[cid2], LDS::EDGE(cid2,cid,tmpD), EC2); // leave?
        if((!intable) && (!no_cache) && (tmpD != LARGER_THAN_UB)){
          insert(cid, cid2, tmpD);
        }
      }
    }
    // return make_pair(-1, -1); 
  }

  // store the closest nn in edges
  // assume edges[cid] already has max written
  // todo: move kdtree to distComputer
  inline void getNN(intT cid, double ub = numeric_limits<double>::max(), intT t_nn = -1){
    // if(edges[cid].getW() ==0){ can't stop, need the one with smallest id
    //   return;
    // }

    // after a round of C <= 50, we might not have all entries
    // in the table. We only have terminal nodes->all clusters
    if(C <= NAIVE_THRESHOLD){
      getNN_naive(cid, ub, t_nn);
      return;
    }

    // if(cid == 23167){
    //   cout << "::: " << t_nn << " " << ub << endl;
    // }

    // check in merge, no inserting merged entry
    // can't writemin to all edges first and then search
    // maybe because a bad neighbor can write to the edge and give a bad radius
    double minD = ub;
    intT nn = t_nn;
    bool intable;

    if(ub == numeric_limits<double>::max()){
        typedef FINDNN::NNsingle<dim, kdnodeT> Fs;
        Fs fs = Fs(uf, cid, eps);
        if(distComputer->nn_process){
          distComputer->getRadius(cid, kdtree->root, &fs);
        }else{
          // treeT treetmp = treeT(PP + cid, 1, false);
          pointT centroid = pointT(getNode(cid)->center, cid);
          treeT treetmp = treeT(&centroid, 1, false);
          // closest to a single point in cluster
          FINDNNP::dualtree<kdnodeT, Fs>(treetmp.root, kdtree->root, &fs); 
        }
        
        nn = uf->find(fs.e->second);
        tie(minD, intable) = getDist(cid, nn); 
        if((!intable) && (!no_cache)  && (minD != LARGER_THAN_UB)){
          insert(cid, nn, minD);
        }
        // if(getNode(nn)->n == 1 && minD == fs.e->getW()){ //useful for early rounds
        //   utils::writeMin(&edges[cid], LDS::EDGE(cid, nn, minD), EC2); 
        //   return ;
        // }
    }

    utils::writeMin(&edges[cid], LDS::EDGE(cid, nn, minD), EC2); 
    utils::writeMin(&edges[nn], LDS::EDGE(nn, cid, minD), EC2);
    // if(cid == 23167){
    //   cout << "::: " << nn << " " << minD << endl;
    // }
    if(minD ==0){
      return;
    }

    Fr fr = Fr(uf, cid, nodes, rootIdx, cacheTbs, edges, distComputer, no_cache, C, eps); 
#ifdef PERF_RANGE
    fr.setCounter(distance_computed, pointsInRange);
#endif  
    // point<dim> pMin1, pMax1; 
    // tie(pMin1, pMax1) = fr.getBox(getNode(cid), minD+eps);
    // FINDNN::NNcandidate<point<dim>, Fr, kdnodeT>(kdtree->root, pMin1, pMax1, &fr);
    double r = fr.getBall(getNode(cid), minD+eps);
    FINDNN::NNcandidate<point<dim>, Fr, kdnodeT>(kdtree->root, getNode(cid)->center, r, &fr);
    // if(cid == 23167){
    //   cout << "done " << fr.getFinalNN() << " "  << fr.getFinalDist() << endl;
    // }

    if(fr.local){
    nn = fr.getFinalNN();
    minD = fr.getFinalDist();
    utils::writeMin(&edges[cid], LDS::EDGE(cid, nn, minD), EC2); 
    }
  }

    //Ave: merge two trees
    // u, v are cluster ids
  inline void merge(intT u, intT v, intT newc, intT round, double height){
    intT rootNodeIdx = nodeIdx.fetch_add(1);
    nodes[rootNodeIdx] = nodeT(newc, round, rootNodeIdx, getNode(u), getNode(v), height);
    rootIdx[newc] = rootNodeIdx;
  }

    // return true only when insert if sucussful
    // neex the  do Swap
    // if one is old and one is new, always check tb[new]->old
    // because the one stored in old might be outdated
    // must have enough space due to alloc new
    // calculate distance if return true
  inline bool insert_check(intT qid, intT rid, bool doSwap){
      CHECK_NO_CACHE(473)
      if(doSwap && qid > rid){
          swap(qid, rid);
      }
      intT qIdx = idx(qid);
      intT rIdx = idx(rid);
      bool inserted; bool reach_thresh;
      tie(inserted, reach_thresh) = cacheTbs[qIdx]->insert_thresh(hashClusterAveET(rid, rIdx, CHECK_TOKEN));
      if(!reach_thresh) return inserted;

      if(doSwap){
        tie(inserted, reach_thresh) = cacheTbs[rIdx]->insert_thresh(hashClusterAveET(qid, qIdx, CHECK_TOKEN));
        if(!reach_thresh) return inserted;
      }
      
      return false;
  }

  inline void updateDist(intT newc){
    CHECK_NO_CACHE(497)
    intT idx1 = leftIdx(newc);
    intT idx2 = rightIdx(newc);
    intT cid1 = getNode(newc)->left->cId;
    intT cid2 = getNode(newc)->right->cId;

    distCacheT *tb1 = cacheTbs[idx1];
    distCacheT *tb2 = cacheTbs[idx2];

    //TODO: optimize to  alloc only once
    _seq<typename distCacheT::eType> TAR1 = tb1->entries();
    _seq<typename distCacheT::eType> TAR2 = tb2->entries();
    // if(C < 100)cout << "sizes in merge: " << TAR1.n+TAR2.n << "/" << tb1->m + tb2->m << endl;
    distCacheT *newtb = cacheTbs[idx(newc)];

    parallel_for(intT i = 0; i<(TAR1.n+TAR2.n); ++i){
      distCacheT::eType *TA = TAR1.A;
      intT offset = 0;
      if(i >= TAR1.n){ //switch to process tb2
        TA = TAR2.A;
        offset = TAR1.n;
      }
      intT j = i-offset;
      intT newc2 = -1;
      intT storedIdx = TA[j].idx; 
      // if storedIdx is inconsistant, the newc2 has changed and we do not calculate
      newc2 = uf->find(TA[j].first);
      if(newc2 != cid1 && newc2 != cid2){

          bool success = false;
          double d;
          // if newc2 is a merged cluster
          if(getNode(newc2)->round == getNode(newc)->round){
            // TA[j].first should == left idx or  right idx
            if(storedIdx == idx(getNode(newc2)->left) || storedIdx == idx(getNode(newc2)->right)){
              success = insert_check(newc, newc2, true); // table might not be symmetric, faster than no ins check
              if(success) d = getNewDistN(newc, newc2); 
            }
            
          }else{
            // TA[j].first should == idx
            if(storedIdx == idx(getNode(newc2))){
              success = insert_check(newc, newc2, false);
              if(success) d = getNewDistO(newc, newc2);
            }
          }
        if(success) { // only insert duplicated entries once 
          insert_helper(newc, newc2, d, newtb, cacheTbs[idx(newc2)]);
        }
      }

    }
    free(TAR1.A);
    free(TAR2.A);
  }

  inline void cleanDist(intT newc){
    CHECK_NO_CACHE(554)
    intT idx1 = leftIdx(newc);
    intT idx2 = rightIdx(newc);
    cacheTbs[idx1]->del();
    cacheTbs[idx2]->del();
    delete cacheTbs[idx1];
    delete cacheTbs[idx2];
    cacheTbs[idx1] = nullptr;
    cacheTbs[idx2] = nullptr;
  }

  inline void checkCache(){
#ifdef BENCHCACHE
      if(no_cache) return;
      parallel_for(intT i = 0; i < C; ++i){ //
        intT cid = activeClusters[i];
          usedCacheCounters[i] = cacheTbs[idx(cid)]->count();
      }
      used_cache=sequence::plusReduce(usedCacheCounters, C);
      long min_cache=sequence::maxIndex(usedCacheCounters, C, less<int>());
      long max_cache=sequence::maxIndex(usedCacheCounters, C, greater<int>());

      cout << "used-cache: " << used_cache << endl;
      cout << "avg-used-cache: " << 1.0*used_cache/C << endl;
      cout << "max-used-cache: " << usedCacheCounters[max_cache] << endl;
      cout << "min-used-cache: " << usedCacheCounters[min_cache] << endl;

      parallel_for(intT i = 0; i < C; ++i){ //
        usedCacheCounters[i] = 0;
      }

#endif   
  }

  // find new activeClusters array based on uf, update C
  inline void updateActiveClusters(intT round){
#ifdef DEBUG
    UTIL::PrintVec(activeClusters, C);
#endif
    
    parallel_for(intT i = 0; i < C; ++i){
      intT cid  = activeClusters[i];
      flag[i] = (uf->find(cid)==cid);
    }

    _seq<intT> newClusters = sequence::pack<intT, intT, sequence::getA<intT,intT> >(NULL, flag, 0, C, sequence::getA<intT,intT>(activeClusters));
    free(activeClusters);
    activeClusters = newClusters.A;
    C = newClusters.n;
    if(marker.doMark(C, round)){ //TODO: move marker to distComputation
      // typedef FINDNN::MarkClusterId<dim, kdnodeT> M;
      // M *cidMarker = new M(uf);
      FINDNN::singletree<kdnodeT, M, typename M::infoT>(kdtree->root, &marker, marker.initVal);
      // delete cidMarker;
    }

#ifdef DEBUG
    cout <<  "pack activeClusters" << endl;
    UTIL::PrintVec(activeClusters, C);
#endif
  }

  void report_perf_range(intT chainNum){
#ifdef PERF_RANGE
  long total_distance_computed =sequence::plusReduce(distance_computed, ELTPERCACHELINE * getWorkers());
  long total_points_in_range=sequence::plusReduce(pointsInRange, ELTPERCACHELINE * getWorkers());

  cout << "distance-computed: " << total_distance_computed << endl;
  cout << "points-in-range: " << total_points_in_range << endl;
  cout << "range-saved: " << (long)chainNum*chainNum -  total_points_in_range << endl;
  cout << "cache-saved: " << total_points_in_range - total_distance_computed << endl;

  for(intT i=0; i<ELTPERCACHELINE * getWorkers(); ++i) {distance_computed[i]=0;pointsInRange[i]=0;}
#endif  
  }

  // edges[i] stores the ith nn  of point i
  // initialize the chain in info
  inline void initChain(TreeChainInfo<dim> *info){
    typedef FINDNN::AllPtsNN<dim, kdnodeT> F;
    F *f = new F(edges, eps);
    FINDNNP::dualtree<kdnodeT, F>(kdtree->root, kdtree->root, f, false);
    parallel_for(intT i=0; i<n; ++i) {
      info->updateChain(edges[i].first, edges[i].second, edges[i].getW());
    }
#ifdef DEBUG
    UTIL::PrintVec2<LDS::EDGE>(edges, n);
#endif
    delete f;
    if(!no_cache){
    parallel_for(intT cid = 0; cid < n; ++cid){
      insert_helper(cid, edges[cid].second,  edges[cid].getW(), cacheTbs[cid], cacheTbs[edges[cid].second]);
    }
    }
  }

// assume table won't be full
inline void fillTable(){
  if(no_cache) return;
  for(intT i=0; i<C; ++i){
    intT cid1 = activeClusters[i];
    for(intT j = i+1; j < C; ++j){
      intT cid2 = activeClusters[j];
      double tmpD; bool intable;
      tie(tmpD, intable) = getDist(cid1, cid2);
      if((!intable)) insert(cid1, cid2, tmpD);
    }
  }
}
  // if((!finder->fullTableFlag) && (finder->C <= 50)){
  //   finder->fillTable();
  //   finder->fullTableFlag = true;
  // }


}; // finder end

namespace ChainTree {
template<int dim, class TF>
inline void chain_find_nn(intT chainNum, TF *finder, TreeChainInfo<dim> *info){
  parallel_for(intT i = 0; i < chainNum; ++i){
      intT cid = info->terminal_nodes[i];
      finder->edges[cid] = LDS::EDGE(-1,-1,numeric_limits<double>::max());
  }

  parallel_for_1(intT i = 0; i < chainNum; ++i){
    intT cid = info->terminal_nodes[i];
    pair<intT, double> prev = info->getChainPrev(i);
    finder->getNN(cid, prev.second, prev.first);
  }
  parallel_for(intT i = 0; i < chainNum; ++i){
      intT cid = info->terminal_nodes[i];
      info->updateChain(cid, finder->edges[cid].second, finder->edges[cid].getW());
#ifdef DEBUG
  UTIL::PrintFunctionItem("NN", "cid", cid);
  UTIL::PrintFunctionItem("NN", "nn", finder->edges[cid].second);
  UTIL::PrintFunctionItem("NN", "w", finder->edges[cid].getW());
#endif
  }
}

template<int dim, class TF>
inline void link_terminal_nodes(UnionFind::ParUF<intT> *uf, TF *finder, TreeChainInfo<dim> *info, intT round, intT *flags){
  intT chainNum = info->chainNum;
  LDS::EDGE *edges = finder->edges;
  
#ifdef TIMING2
  timer t1; t1.start();
#endif

  parallel_for(intT i = 0; i < chainNum; ++i){
    intT cid = info->terminal_nodes[i];
    // intT cid = uf->find(edges[i].first); 
    intT nn = uf->find(edges[cid].second);
    intT nn_of_nn =info->getNN(nn);
    // if( nn_of_nn> 0) nn_of_nn = uf->find(nn_of_nn); do not need because all clusters on chain did not merge in the last round

    // avoid merging twice
    // update chain info and merge trees
    if ((!(info->findNN[nn] && cid <= nn)) && (nn != cid && nn_of_nn == cid)){// if they are RNN
      intT newc = uf->link(cid, nn, edges[cid].getW());
      if(newc == cid){
        info->invalidate(nn, -2);
        info->invalidate(cid, -1);
      }else{
        info->invalidate(cid, -2);
        info->invalidate(nn, -1);
      }
      info->invalidateRev(newc);
      finder->merge(cid, nn, newc, round, edges[cid].getW());
      flags[i] = newc;
    }else{
      flags[i] = -1;
    }
  }
  
  if(!finder->no_cache){
  _seq<intT> merged = sequence::filter(flags, chainNum, [&](intT i){return i!=-1;});

#ifdef TIMING2
	 if(LINKAGE_DOPRINT(round)){ UTIL::PrintSubtimer(":::merge", t1.next());  cout << "merged.n: " << merged.n << endl;}
#endif
  // insert to new hashtables and delete old hashtables
  parallel_for(intT i = 0; i < merged.n; ++i){
    intT newc = merged.A[i];//flags[i];
    // if(newc != -1){
      finder->updateDist(newc);
    // }
  }
#ifdef TIMING2
	if(LINKAGE_DOPRINT(round)){ UTIL::PrintSubtimer(":::update dist", t1.next());}
#endif

#ifndef ALLOCDOUBLE
  // delete old hashtables
  parallel_for(intT i = 0; i < merged.n; ++i){
    intT newc = merged.A[i];//flags[i];
    // if(newc != -1){
      finder->cleanDist(newc);
    // }
  }
#ifdef TIMING2
	if(LINKAGE_DOPRINT(round)){ UTIL::PrintSubtimer(":::delete old table", t1.next());}
#endif
#endif
  free(merged.A);
  } //  end !no_cache

}

// TODO: proceed in stages to save work
template<int dim, class TF>
inline void chain_linkage(TF *finder, timer t1){
  UTIL::PrintCaption("CHAIN TREE");
  cout << "num workers " << getWorkers() << endl;
  // cout << "MAX_CACHE_TABLE_SIZE " << MAX_CACHE_TABLE_SIZE << endl;
  cout << "hash table size " << finder->hashTableSize << endl;
  
  UnionFind::ParUF<intT> *uf = finder->uf;
  intT n = finder->n;
  intT chainNum = n;
  TreeChainInfo<dim> *info = new TreeChainInfo<dim>(n, finder->eps);
  intT *flags = newA(intT, chainNum);// used in linkterminalnodes

#ifdef VERBOSE
	UTIL::PrintSubtimer("initialize", t1.next());
#endif
//  ofstream file_obj;
//  file_obj.open("debug/avg_UCI1_32_1th.txt"); //"+ to_string(round) + "

  int round = 0;
  bool print = false;
  while(finder->C > 1 ){
    round ++;
#ifdef VERBOSE
    if(LINKAGE_DOPRINT(round)){//
    print = true;
    UTIL::PrintBreak();
    UTIL::PrintFunctionItem("CLINK", "Round", round);
    UTIL::PrintFunctionItem("CLINK", "Comp Num", finder->C);
    UTIL::PrintFunctionItem("Chain", "#", info->chainNum);}else{print = false;}
#endif
  if(round >= n)  exit(1);
  if(round >= 2 && info->chainNum == 0) zero_chain_debug(finder, round, info); 

  if(round == 1){
    finder->initChain(info);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("init-chains", t1.next());
#endif
  }else{
    chain_find_nn(chainNum, finder, info);
    // findAllNNBruteForce(chainNum, finder, info);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("find-nn", t1.next());
#endif
  }

#ifdef BENCHCACHE
      finder->checkCache();
      if(print) UTIL::PrintSubtimer("check-cache", t1.next());
#endif
#ifdef PERF_RANGE
    finder->report_perf_range(info->chainNum);
    t1.next();
#endif
  // cout << " ----" << endl;
  // cout << finder->edges[41705].first << " " << finder->edges[41705].second << " "  << finder->edges[41705].getW() << endl;
//  UTIL::printTuple(finder->getDist(11791, 1470));
//  cout << finder->getDistNaiveDebug(11791, 1470) << endl;
  // cout << " ----" << endl;
//  UTIL::PrintTable(finder->getTable(11791));

// if(true){
//    for(intT i = 0; i < finder->C; ++i){
//      file_obj << round << " " << finder->activeClusters[i] << endl;
//    }
//    file_obj << round << "========" << endl;
//     for(intT i = 0; i < chainNum; ++i){
//       intT cid = info->terminal_nodes[i];
//       file_obj << round << " " << cid << " " << finder->edges[cid].second << " " << finder->edges[cid].getW() << endl;
//     }

//    file_obj << round << "========" << endl;
// } //end true

    link_terminal_nodes<dim, TF>(uf, finder, info, round, flags);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("link-update", t1.next());
#endif
      // get ready for next round
    finder->updateActiveClusters(round);
    finder->distComputer->update(round, finder);
    info->next(finder);
    chainNum = info->chainNum;

#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("update-clusters", t1.next());
#endif
  t1.next();
  }
  UTIL::PrintFunctionItem("CLINK", "rounds", round);
  delete info;
  free(flags);
  finder->distComputer->postProcess(finder);
}

}//end of namespace
#endif
