#ifndef METHOD_BORUVKA_HDBSCAN_H
#define METHOD_BORUVKA_HDBSCAN_H

#include "sequence.h"
#include "geometry.h"
#include "sampleSort.h"
#include "gettime.h"
#include "unionfind.h"

#include "kdTree2.h"
#include "kNearestNeighbors.h"
#include "shared.h"
#include "dist.h"
#include "neighbor.h"
#include "neighbor_parallel.h"

#define PRINT_FREQ 1
#define LINKAGE_DOPRINT(round) (round < 5 || round % PRINT_FREQ == 0)

#include <limits>
#include <cstdint>

using namespace std;

//TODO: add nodoeIdx to Info 
template<int dim, class distF,  class nodeInfo>
class BoruvkaFinder { //: public NNFinder<dim>
  public:
  typedef iPoint<dim> pointT;
  // typedef typename FINDNN::CLinkNodeInfo nodeInfo;
  typedef FINDNN::LinkNodeInfo1<dim, pointT *> nodeT;
  typedef kdNode<dim, pointT, nodeInfo> kdnodeT;
  typedef kdTree<dim, pointT, nodeInfo > treeT;
  typedef FINDNN::AllPtsNNMetric<dim, kdnodeT, distF> FNN;
  typedef FINDNN::MarkClusterIdResetUpperDiam<dim, FINDNN::DummyRangeF<dim, pointT, nodeInfo, nodeT>, distF> M;

//   bool no_cache;
  intT C;// number of clusters
  intT n;// number of points
  pointT *PP;
  intT *rootIdx;
  UnionFind::ParUF<intT> *uf;
  intT *activeClusters; // ids of connected components
  intT *activeClusters2; // ids of connected components, used to swap
//   intT *CMap; // map from cluster id to location in activeClusters;
  intT *flags;//used in updateActivateClusters
  treeT *kdtree;
  LDS::EDGE *edges; //edges[i] stores the min neigh of cluster i

  nodeT *nodes; // preallocated space to store tree nodes
  const LDS::edgeComparator2 EC2 = LDS::edgeComparator2();
  atomic<intT> nodeIdx; // the index of next node to use for cluster trees

  distF *distComputer;
  FNN fnn;
  M marker;

  BoruvkaFinder(intT t_n, point<dim>* t_P, UnionFind::ParUF<intT> *t_uf, commandLine params): uf(t_uf), n(t_n){
    C = n;
    activeClusters = newA(intT, 2*n);
    activeClusters2 = activeClusters+n;
    parallel_for(intT i=0; i<n; ++i) {activeClusters[i] = i;}
    // CMap = newA(intT, n);
    // parallel_for(intT i=0; i<n; ++i) {CMap[i] = i;}
    PP = newA(pointT, n);
    parallel_for(intT i=0; i<n; ++i) {PP[i] = pointT(t_P[i], i);}

    distComputer = new distF(PP, n, params);

    rootIdx = newA(intT, n);
    parallel_for(intT i=0; i<n; ++i) {rootIdx[i] = i;}
    nodes = newA(nodeT, 2*n);
    // distComputer->initNodes(nodes, n);
    
    flags = newA(intT, C);

    edges = (LDS::EDGE *)malloc(sizeof(LDS::EDGE) * n);
    parallel_for(intT i=0; i<n; ++i) {edges[i] = LDS::EDGE(-1,-1,numeric_limits<double>::max());}

    nodeIdx.store(n); // have used n nodes
    kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);
    distComputer->template initTree<treeT>(kdtree, n);
        
    fnn = FNN(edges, uf, distComputer);
    marker = M(uf, nodes, rootIdx, distComputer);
  }

  ~BoruvkaFinder(){
    free(PP);
    free(flags);
    if(reinterpret_cast<std::uintptr_t>(activeClusters) < reinterpret_cast<std::uintptr_t>(activeClusters2)){
      free(activeClusters);
    }else{
      free(activeClusters2);
    }
    // free(CMap);
    delete kdtree;
    free(nodes);
    free(edges);
  }

  inline intT cid(nodeT* node){ return node->cId;}
  inline intT idx(nodeT* node){ return node->idx;}
  inline nodeT *getNode(intT cid){return nodes+rootIdx[cid];}
  inline intT idx(intT cid){return idx(getNode(cid));}
  inline intT leftIdx(intT cid){return getNode(cid)->left->idx;}
  inline intT rightIdx(intT cid){return getNode(cid)->right->idx;}
  inline intT cid(intT idx){return cid(nodes[idx]);}

  struct ClusterIdComparator {
  UnionFind::ParUF<intT> *uf;
  ClusterIdComparator(UnionFind::ParUF<intT> *t_uf):uf(t_uf) {}

  bool operator() (const intT &i, const intT &j) {
    return uf->find(i) < uf->find(j);
  }
  };

  inline void link_terminal_nodes(intT round){

    parallel_for(intT i = 0; i < C; ++i){
        intT cid = activeClusters[i];
        intT nn = edges[cid].second;//uf->find();
        uf->link(cid, nn, edges[cid].getW());
    }

    sampleSort(activeClusters, C, ClusterIdComparator(uf));
    flags[0] = 1;
    parallel_for (intT i = 1; i < C; ++ i) {
        if (uf->find(activeClusters[i]) != uf->find(activeClusters[i-1])) {
        flags[i] = 1;
        }else{
        flags[i] = 0;
        }
    }
    intT numClusters;
    numClusters = sequence::prefixSum<intT>(flags, 0, C);
    intT *clusterOffsets = newA(intT, numClusters + 1);
    parallel_for (intT i = 0; i < C - 1; ++ i) {
        if (flags[i] != flags[i + 1]) {
        clusterOffsets[flags[i]] = i;
        }
    }
    intT i = C - 1;
    if (flags[i] != numClusters) {
        clusterOffsets[flags[i]] = i;
    }
    clusterOffsets[numClusters] = C;

    parallel_for (intT i = 0; i < numClusters; ++ i) {
        activeClusters2[i] = uf->find(activeClusters[clusterOffsets[i]]);
    }

    /////// maintaining hierarchy
    intT nodeIdxBase = nodeIdx.load()-1;
    parallel_for (intT i = 0; i < numClusters; ++ i) {
        intT newc = activeClusters2[i];
        parallel_for (intT j = clusterOffsets[i]+1; j < clusterOffsets[i+1]; ++ j) {
            intT newNodeIdx = nodeIdxBase + j - i;
            nodeT *left = nodes+newNodeIdx-1;
            if(j == clusterOffsets[i]+1){
                left = getNode(activeClusters[j-1]);
            }
            nodes[newNodeIdx] = nodeT(newc, round, newNodeIdx, left, getNode(activeClusters[j]));
        }
        rootIdx[newc] = nodeIdxBase + clusterOffsets[i+1] - i;
    }

        nodeIdx.fetch_add(C-numClusters);
    /////// maintaining hierarchy

        swap(activeClusters, activeClusters2);
        C = numClusters;
        // parallel_for(intT i = 0; i < C; ++i){
        //   intT cid  = activeClusters[i];
        //   CMap[cid] = i;
        // }
        if(marker.doMark(C, round)){ //TODO: move marker to distComputation
        FINDNN::singletree<kdnodeT, M, typename M::infoT>(kdtree->root, &marker, marker.initVal);
        }
    }

    // init kdtree upperbounds using diameter
    // void init_upper_bound(kdnodeT *Q) {
    //   // Q->nInfo.resetUB(); need to add id not done in marking
    //   if (!Q->isLeaf()) {
    //     if (Q->size() > 2000) {
    //       cilk_spawn init_upper_bound(Q->L());
    //       init_upper_bound(Q->R());
    //       cilk_sync;
    //     } else {
    //       init_upper_bound(Q->L());
    //       init_upper_bound(Q->R());
    //     }
    //     (Q->nInfo).updateUB(max((Q->left->nInfo).getUB(), (Q->right->nInfo).getUB()));
    //   }

    //   // an optimization for AllOneBCP, where 1-BCP edge is computed for all clusters
    //   // when no known upperbound is known, use diameter of node
    //   // only applicable when multiple clusters are in the node (m_cluster_id)
    //   if ( Q->nInfo.getCId() < 0){
    //     (Q->nInfo).updateUB(Q->Diameter());
    //   }
    // }

  inline void findAllNN(intT round){
      if(distComputer->useSimpleNN(round)){ //using a simpler findNN
        FINDNN::AllPtsNN<dim, kdnodeT> *f = new FINDNN::AllPtsNN<dim, kdnodeT>(edges);
        FINDNNP::dualtree<kdnodeT, FINDNN::AllPtsNN<dim, kdnodeT>>(kdtree->root, kdtree->root, f, false);
        delete f;
        return;
      }
    parallel_for(intT i = 0; i < C; ++i){
      intT cid = activeClusters[i];
      edges[cid] = LDS::EDGE(-1,-1,numeric_limits<double>::max());
    }
    FINDNNP::dualtree<kdnodeT, FNN>(kdtree->root, kdtree->root, &fnn, false);
  }

}; // finder end


// TODO: proceed in stages to save work
template<int dim, class TF>
inline void boruvka_linkage(TF *finder, timer t1){
  UTIL::PrintCaption("Boruvka Method");
  cout << "num workers " << getWorkers() << endl;
  
  UnionFind::ParUF<intT> *uf = finder->uf;
  intT n = finder->n;

#ifdef VERBOSE
	UTIL::PrintSubtimer("initialize", t1.next());
#endif
//  ofstream file_obj;
//  file_obj.open("debug/single_19K_wrong.txt"); //"+ to_string(round) + "

  int round = 0;
  bool print = false;
  while(finder->C > 1 ){
    round ++;
#ifdef VERBOSE
    if(LINKAGE_DOPRINT(round)){//
    print = true;
    UTIL::PrintBreak();
    UTIL::PrintFunctionItem("CLINK", "Round", round);
    UTIL::PrintFunctionItem("CLINK", "Comp Num", finder->C);}else{print = false;}
#endif
    // TODO: add k-nn finding step?
  if(round >= n)  exit(1);
    finder->findAllNN(round);
    // finder->findAllNNBruteForce();
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("find-nn", t1.next());
#endif
    finder->link_terminal_nodes(round);
    // finder->distComputer->update(round, finder);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("update-clusters", t1.next());
#endif
  t1.next();
  }
  UTIL::PrintFunctionItem("CLINK", "rounds", round);

}
#endif


// if(round > 1){
//    for(intT i = 0; i < finder->C; ++i){
//      file_obj << round << " " << finder->activeClusters[i] << endl;
//    }
//    file_obj << round << "========" << endl;
//     for(intT i = 0; i < finder->C; ++i){
//       intT cid = finder->activeClusters[i];
//       file_obj << round << " " << cid << " " << finder->edges[cid].second << " " << finder->edges[cid].getW() << endl;
//       //timer t2;t2.start();
//       // double dist_brute = finder->distComputer->getDistNaive(finder->getNode(cid), finder->getNode(finder->edges[cid].second));
//       //cout << t2.next() << endl;
//       // if(abs(finder->edges[cid].getW()  - dist_brute) > 0.000001){
//       //     cout << std::fixed;
//       //     cout << std::setprecision(10);
//       //     cout << "exit " << round << " " << cid << " " << finder->edges[cid].second <<" "<< dist_brute << " " << finder->edges[cid].getW() << endl;
//  	    //     cout << finder->find(cid, finder->edges[cid].second) << endl;
//  	    //     exit(1);
//       // }
//     }

//    file_obj << round << "========" << endl;
// } //end round>1



  // inline void findAllNNBruteForce(){
  //       parallel_for(intT i = 0; i < C; ++i){
  //           intT cid = activeClusters[i];
  //           edges[cid] = LDS::EDGE(-1,-1,numeric_limits<double>::max());
  //       }
  //       for(intT i = 0; i< n; ++i){
  //           for(intT j = i+1; j< n; ++j){
  //               intT cid1 = uf->find(i);
  //               intT cid2 = uf->find(j);
  //               if(cid1 == cid2) continue;
  //               double qrdist = PP[i].pointDist(PP[j]);
  //               qrdist = max(qrdist, max(distComputer->coredists[PP[i].i], distComputer->coredists[PP[j].i]));
  //               utils::writeMin(&edges[cid1], LDS::EDGE(cid1,cid2,qrdist), EC2);
  //               utils::writeMin(&edges[cid2], LDS::EDGE(cid2,cid1,qrdist), EC2);
  //           }
  //       }
  //   }