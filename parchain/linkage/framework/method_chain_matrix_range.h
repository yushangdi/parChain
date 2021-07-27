#ifndef METHOD_CHAIN_MATRIX_RANGE_H
#define METHOD_CHAIN_MATRIX_RANGE_H

#include "sequence.h"
#include "geometry.h"
#include "sampleSort.h"
#include "gettime.h"
#include "unionfind.h"

#include "kdTree2.h"

#include <limits>

#include "shared.h"
#include "chain.h"
#include "cluster.h"
#include "dist.h"
#include "neighbor.h"
#include "neighbor_parallel.h"
#include "rangeQuery.h"

#include "distMatrix.h"
#include "dendrogram.h"
using namespace std;

//TODO: add nodoeIdx to Info 
template<int dim, class distF, class Fr>
class MatrixRangeNNFinder { //: public NNFinder<dim>
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

  intT C;// number of clusters
  intT n;// number of points
  pointT *PP;
  intT *rootIdx;
  UnionFind::ParUF<intT> *uf;
  intT *activeClusters; // ids of connected components
  bool *flag;//used in updateActivateClusters
  treeT *kdtree;
  LDS::EDGE *edges; //edges[i] stores the min neigh of cluster i

  DM<dim> *matrix; // distance to clusters, store two copies of each distance

  nodeT *nodes; // preallocated space to store tree nodes
  LDS::edgeComparator2 EC2;
  double eps;
  atomic<intT> nodeIdx; // the index of next node to use for cluster trees

  distF *distComputer;
  // M marker;

  intT NAIVE_THRESHOLD = 50;


#ifdef PERF_RANGE
  long *distance_computed;
  long *pointsInRange;
  // long max_max_cache_used_allrounds=0;
  // long max_avg_cache_used_allrounds = 0;
  // long times_cache_full;
#endif  

  //TODO: batch allocate cache table space
  //TODO: free cache table when merges
  MatrixRangeNNFinder(intT t_n, point<dim>* t_P, UnionFind::ParUF<intT> *t_uf, bool dummy_no_cache, double t_eps, intT t_naive_thresh = 5): 
    uf(t_uf), n(t_n), eps(t_eps), NAIVE_THRESHOLD(t_naive_thresh){
    EC2 = LDS::edgeComparator2(eps);
    C = n;
    activeClusters = newA(intT, n);
    parallel_for(intT i=0; i<n; ++i) {activeClusters[i] = i;}
    PP = newA(pointT, n);
    parallel_for(intT i=0; i<n; ++i) {PP[i] = pointT(t_P[i], i);}

    distComputer = new distF(PP, n, /* use_range */ true);

    rootIdx = newA(intT, n);
    parallel_for(intT i=0; i<n; ++i) {rootIdx[i] = i;}
    nodes = newA(nodeT, 2*n);
    distComputer->initNodes(nodes, n);
    
    flag = newA(bool, C);

    edges = (LDS::EDGE *)malloc(sizeof(LDS::EDGE) * n);
    parallel_for(intT i=0; i<n; ++i) {edges[i] = LDS::EDGE(-1,-1,numeric_limits<double>::max());}


    if(distComputer->method == LDS::AVGSQ){
      matrix = new DM<dim>(t_P, n, true);
    }else{
      matrix = new DM<dim>(t_P, n);
    }

    nodeIdx.store(n); // have used n nodes
    kdtree = new treeT(PP, n);
    // marker = M(uf, nodes, rootIdx);

  }

  ~MatrixRangeNNFinder(){
    free(PP);
    free(flag);
    free(activeClusters);
    delete kdtree;
    free(nodes);
    free(edges);
    free(rootIdx);
    delete matrix;
    delete distComputer;
  }

  inline intT cid(nodeT* node){ return node->cId;}
  inline intT idx(nodeT* node){ return node->idx;}
  inline nodeT *getNode(intT cid){return nodes+rootIdx[cid];}
  inline intT idx(intT cid){return idx(getNode(cid));}
  inline intT leftIdx(intT cid){return getNode(cid)->left->idx;}
  inline intT rightIdx(intT cid){return getNode(cid)->right->idx;}
  inline intT cid(intT idx){return cid(nodes[idx]);}

  // i, j are cluster ids
  // find distance in matrix
  tuple<double, bool> getDist(intT i,  intT j, double lb = -1, double ub = numeric_limits<double>::max(), bool par = true){
    return make_tuple(matrix->get(i, j), true);
  }

  tuple<double, bool> getDist(nodeT *inode,  nodeT *jnode, double lb = -1, double ub = numeric_limits<double>::max(), bool par = true){
    return getDist(cid(inode), cid(jnode));
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

    double d1,d2; bool intable;
    tie(d1, intable) = getDist(ql,rroot);

    tie(d2, intable) = getDist(qr,rroot);

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
     intT nn = sequence::maxIndex<pair<double, intT>>(0,C, neighborComparator(eps),neighborDistGetter<DM<dim>>(cid, matrix, activeClusters));
    nn = activeClusters[nn];
    double minD = matrix->get(cid, nn);
    edges[cid].update(cid, nn, minD);
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

    // check in merge, no inserting merged entry
    // can't writemin to all edges first and then search
    // maybe because a bad neighbor can write to the edge and give a bad radius
    double minD = ub;
    intT nn = t_nn;
    bool intable;

    if(ub == numeric_limits<double>::max()){
        typedef FINDNN::NNsingle<dim, kdnodeT> Fs;
        Fs fs = Fs(uf, cid, eps);

          pointT centroid = pointT(getNode(cid)->center, cid);
          treeT treetmp = treeT(&centroid, 1, false);
          // closest to a single point in cluster
          FINDNNP::dualtree<kdnodeT, Fs>(treetmp.root, kdtree->root, &fs); 
        
        nn = uf->find(fs.e->second);
        tie(minD, intable) = getDist(cid, nn); 
    }

    utils::writeMin(&edges[cid], LDS::EDGE(cid, nn, minD), EC2); 
    utils::writeMin(&edges[nn], LDS::EDGE(nn, cid, minD), EC2);
    // if(cid == 23167){
    //   cout << "::: " << nn << " " << minD << endl;
    // }
    if(minD ==0){
      return;
    }

    Fr fr = Fr(uf, cid, nodes, rootIdx, matrix, edges, distComputer, true, C, eps); 
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



  inline void updateDist(intT newc, intT round){
    intT idx1 = leftIdx(newc);
    intT idx2 = rightIdx(newc);
    intT cid1 = getNode(newc)->left->cId;
    intT cid2 = getNode(newc)->right->cId;

    intT oldc = newc == cid1? cid2:cid1;
    parallel_for(intT i = 0; i<C; ++i){
        intT updc = activeClusters[i];
        if(updc != cid1 && updc != cid2 && uf->find(updc) == updc){
          if(getNode(updc)->round == round){ /* if updc is also a merged cluster */
            if(updc > newc){ /* the smaller id cluster will update */
                double dist = getNewDistN(newc, updc); 
                matrix->update(updc, newc, dist);
            }
          }else{
              double dist = getNewDistO(newc, updc);
              matrix->update(updc, newc, dist);
          }

        }//end of (updc != cid1 && updc != cid2)
    }//end parallel for
  }

  // find new activeClusters array based on uf, update C
  inline void updateActiveClusters(intT round){
    parallel_for(intT i = 0; i < C; ++i){
      intT cid  = activeClusters[i];
      flag[i] = (uf->find(cid)==cid);
    }

    _seq<intT> newClusters = sequence::pack<intT, intT, sequence::getA<intT,intT> >(NULL, flag, 0, C, sequence::getA<intT,intT>(activeClusters));
    free(activeClusters);
    activeClusters = newClusters.A;
    C = newClusters.n;
    // if(marker.doMark(C, round)){ 
    //   FINDNN::singletree<kdnodeT, M, typename M::infoT>(kdtree->root, &marker, marker.initVal);
    // }
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
  }

}; // finder end

namespace ChainMatrixRange {
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
  
  _seq<intT> merged = sequence::filter(flags, chainNum, [&](intT i){return i!=-1;});
  // insert to new hashtables and delete old hashtables
  parallel_for(intT i = 0; i < merged.n; ++i){
    intT newc = merged.A[i];
    finder->updateDist(newc, round);
  }
#ifdef TIMING2
	if(LINKAGE_DOPRINT(round)){ UTIL::PrintSubtimer(":::delete old table", t1.next());}
#endif
  free(merged.A);

}

// TODO: proceed in stages to save work
template<int dim, class TF>
inline void chain_linkage(TF *finder, timer t1){
  UTIL::PrintCaption("CHAIN Matrix Range");
  cout << "num workers " << getWorkers() << endl;
  
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

template<int dim, class pointT>
inline void completeLinkage(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps, intT naive_thresh){
  timer t1;
  t1.start();
    // DISTANCE MATRIX
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;
  using nodeInfo = FINDNN::CLinkNodeInfo; // for kdtree

  cout << "Matrix Range Distance Complete Linkage of " << n << ", dim " << dim << " points" << endl; 
  using distT = MatrixDistanceComputer::distComplete3<dim, pointT, nodeT, nodeInfo>;
  using boxT = FINDNN::queryBallSimple<dim, nodeT>; 
  using Fr = FINDNN::RangeQueryMatrixCenterF<dim, pointT, nodeInfo, distT, boxT>;

  using FinderT = MatrixRangeNNFinder<dim, distT, Fr>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, false, eps, naive_thresh); 
  chain_linkage<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());

}

template<int dim, class pointT>
inline void wardLinkage(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps, intT naive_thresh){
  timer t1;
  t1.start();
    // DISTANCE MATRIX
  using nodeT = FINDNN::LinkNodeInfo3<dim, pointT *>;
  using nodeInfo = FINDNN::WLinkNodeInfo; // for kdtree

  cout << "Matrix Range Distance Ward Linkage of " << n << ", dim " << dim << " points" << endl;
  using MCenter = FINDNN::MarkKdTreeCenters<dim, pointT, nodeInfo>;
  using distT = MatrixDistanceComputer::distWard1<dim, pointT, nodeT, nodeInfo, MCenter>;
  using boxT = FINDNN::queryBallWard<dim, nodeT>;
  using Fr = FINDNN::RangeQueryMatrixCenterF<dim, pointT, nodeInfo, distT, boxT>;

  using FinderT = MatrixRangeNNFinder<dim, distT, Fr>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, false, eps, naive_thresh); 
  chain_linkage<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());

}

template<int dim, class pointT>
inline void avgLinkage(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps, intT naive_thresh){
  timer t1;
  t1.start();
    // DISTANCE MATRIX
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;
  using nodeInfo = FINDNN::CLinkNodeInfo; // for kdtree

  cout << "Matrix Range Distance Average Square Linkage of " << n << ", dim " << dim << " points" << endl; 
  // using distT = distAverage1<dim, pointT, nodeT>;
  using distT = MatrixDistanceComputer::distAverage5<dim, pointT, nodeT>;
  using boxT = FINDNN::queryBallSimple<dim, nodeT>; 
  using Fr = FINDNN::RangeQueryMatrixCenterF<dim, pointT, nodeInfo, distT, boxT>;

  using FinderT = MatrixRangeNNFinder<dim, distT, Fr>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, false, eps, naive_thresh); 
  chain_linkage<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());

}

template<int dim, class pointT>
inline void avgsqLinkage(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps, intT naive_thresh){
  timer t1;
  t1.start();
    // DISTANCE MATRIX
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;
  using nodeInfo = FINDNN::CLinkNodeInfo; // for kdtree

  cout << "Matrix Range Distance Average Square Linkage of " << n << ", dim " << dim << " points" << endl; 
  using nodeInfo = FINDNN::CLinkNodeInfo; // for kdtree
  using distT = MatrixDistanceComputer::distAverage4<dim, pointT, nodeT, nodeInfo>; // consecutive point array
  using boxT = FINDNN::queryBallSimple<dim, nodeT>; //still use simple, distT will  postprocess
  using Fr = FINDNN::RangeQueryMatrixCenterF<dim, pointT, nodeInfo, distT, boxT>;

  using FinderT = MatrixRangeNNFinder<dim, distT, Fr>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, false, eps, naive_thresh); 
  chain_linkage<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());

}

}//end of namespace
#endif
