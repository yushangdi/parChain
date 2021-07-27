#ifndef METHOD_CHAIN_MATRIX_H
#define METHOD_CHAIN_MATRIX_H

#include "sequence.h"
#include "geometry.h"
#include "sampleSort.h"
#include "gettime.h"
#include "unionfind.h"

// #define A_HASH_LINKAGE_PROB_PRINT

#include "ndHash.h"
#include "serialHash.h"

#include "kdTree2.h"
#include "shared.h"
#include "chain.h"
#include "cluster.h"
#include "dist.h"
#include "neighbor.h"
#include "neighbor_parallel.h"
#include "rangeQuery.h"

#include <limits>
#include "distMatrix.h"

using namespace std;

//TODO: add nodoeIdx to Info 
template<int dim, class distF>
class MatrixNNFinder {
  public:
  typedef typename distF::pointT pointT;
  typedef typename distF::nodeT nodeT;

  bool no_cache;
  intT C;// number of clusters
  intT n;// number of points
  pointT *PP;
  intT *rootIdx;
  UnionFind::ParUF<intT> *uf;
  intT *activeClusters; // ids of connected components
//   bool *activeFlags; // activeFlags[i] is cluster id i
  bool *flag;//used in updateActivateClusters, flag[i] is activeClusters[i]
  LDS::EDGE *edges; //edges[i] stores the min neigh of cluster i

  DM<dim> *matrix; // distance to clusters, store two copies of each distance
  nodeT *nodes; // preallocated space to store tree nodes
  LDS::edgeComparator2 EC2;
  double eps;
  atomic<intT> nodeIdx; // the index of next node to use for cluster trees

  distF *distComputer;

  
  //TODO: batch allocate cache table space
  //TODO: free cache table when merges
  MatrixNNFinder(intT t_n, point<dim>* t_P, UnionFind::ParUF<intT> *t_uf, 
    bool t_noCache, double t_eps): uf(t_uf), n(t_n), eps(t_eps){
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
    // activeFlags = newA(bool, C);
    // parallel_for(intT i=0; i<C; ++i) {activeFlags[i] = true;}

    edges = (LDS::EDGE *)malloc(sizeof(LDS::EDGE) * n);
    parallel_for(intT i=0; i<n; ++i) {edges[i] = LDS::EDGE(-1,-1,numeric_limits<double>::max());}

    nodeIdx.store(n); // have used n nodes

    if(distComputer->method == LDS::AVGSQ){
      matrix = new DM<dim>(t_P, n, true);
    }else{
      matrix = new DM<dim>(t_P, n);
    }
    

  }

  ~MatrixNNFinder(){
    free(PP);
    free(flag);
    // free(activeFlags);
    free(activeClusters);
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
    double dij = getNode(newc)->getHeight();
    nodeT* rroot = getNode(rid);
    intT nr = rroot->n;

    double d1,d2; bool intable;
    tie(d1, intable) = getDist(ql,rroot);

    tie(d2, intable) = getDist(qr,rroot);

    return distComputer->updateDistO(d1, d2, nql, nqr, nr, dij);
  }


  //newc is a newly merged cluster
  // newc is new, rid is merged
  //todo: cids on cluster tree do not need to be marked
  // overwrite entry?
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
  inline void getNN(intT cid){
    intT nn = sequence::maxIndex<pair<double, intT>>(0,C, neighborComparator(eps),neighborDistGetter<DM<dim>>(cid, matrix, activeClusters));
    nn = activeClusters[nn];
    double minD = matrix->get(cid, nn);
    edges[cid].update(cid, nn, minD);
    // utils::writeMin(&edges[cid], LDS::EDGE(cid, nn, minD), EC2); 
    // utils::writeMin(&edges[nn], LDS::EDGE(nn, cid, minD), EC2);

  }

    //Ave: merge two trees
    // u, v are cluster ids
  inline void merge(intT u, intT v, intT newc, intT round, double height){
    intT rootNodeIdx = nodeIdx.fetch_add(1);
    nodes[rootNodeIdx] = nodeT(newc, round, rootNodeIdx, getNode(u), getNode(v), height);
    rootIdx[newc] = rootNodeIdx;
    // if(newc==u){
    //     activeFlags[v] = false;
    // }else{
    //     activeFlags[u] = false;
    // }
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
  }
}; // finder end



namespace ChainMatrix {
template<int dim, class TF>
inline void chain_find_nn(intT chainNum, TF *finder, TreeChainInfo<dim> *info){
    // update edges array
  parallel_for_1(intT i = 0; i < chainNum; ++i){
    intT cid = info->terminal_nodes[i];
    finder->getNN(cid);
  }
  // update chain
  parallel_for(intT i = 0; i < chainNum; ++i){
      intT cid = info->terminal_nodes[i];
      info->updateChain(cid, finder->edges[cid].second, finder->edges[cid].getW());
  }
}

//TODO: cache idx and cids
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

#ifdef TIMING2
	 if(LINKAGE_DOPRINT(round)){ UTIL::PrintSubtimer(":::merge", t1.next());  cout << "merged.n: " << merged.n << endl;}
#endif

  // insert to new hashtables and delete old hashtables
  parallel_for(intT i = 0; i < merged.n; ++i){
    intT newc = merged.A[i];
    finder->updateDist(newc, round);
  }
#ifdef TIMING2
	if(LINKAGE_DOPRINT(round)){ UTIL::PrintSubtimer(":::update dist", t1.next());}
#endif
  free(merged.A);

}

// TODO: proceed in stages to save work
template<int dim, class TF>
inline void chain_linkage_matrix(TF *finder, timer t1){
  UTIL::PrintCaption("CHAIN MATRIX");
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
//  file_obj.open("debug/avg_2DVar1M_edges_96.txt"); //"+ to_string(round) + "

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
    chain_find_nn(chainNum, finder, info);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("find-nn", t1.next());
#endif
    link_terminal_nodes<dim, TF>(uf, finder, info, round, flags);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("link-update", t1.next());
#endif
      // get ready for next round
    finder->updateActiveClusters(round);
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
  // finder->distComputer->postProcess(finder);
}

}//end pf namespace
#endif