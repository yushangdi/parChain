#ifndef METHOD_NAIVE_TREE_H
#define METHOD_NAIVE_TREE_H

#include "shared.h"
#include "kdTree2.h"
#include "method_chain_tree.h"

using namespace std;

// *************************************************************
//   Naive Linkage with Tree
// *************************************************************



template<int dim>
inline void naive_find_nn(LDS::EDGE *edges, TreeNNFinder<dim> *finder){

  parallel_for(intT i = 0; i < finder->C; ++i){
    intT cid = finder->activeClusters[i];
    pair<intT, double> nnDist = finder->getNN_naive(cid);//TODO: add lower bound = diameter
    edges[i] = LDS::EDGE(cid, nnDist.first, nnDist.second);
// #ifdef DEBUG
//     UTIL::PrintFunctionItem("NN", "cid", cid);
//     UTIL::PrintFunctionItem("NN", "nn", nnDist.first);
// #endif
    }
}

template<int dim>
inline void naive_find_nn_range(LDS::EDGE *edges, TreeNNFinder<dim> *finder){
  // intT PNum =  getWorkers(); 
  // intT *sizes = newA(intT, PNum * finder->n);//TODO: optimize to C size

  parallel_for(intT i = 0; i < finder->C; ++i){
    intT cid = finder->activeClusters[i];
    pair<intT, double> nnDist = finder->getNN(cid); //, numeric_limits<double>::max(), -1
    edges[i] = LDS::EDGE(cid, nnDist.first, nnDist.second);
// #ifdef DEBUG
//     UTIL::PrintFunctionItem("NN", "cid", cid);
//     UTIL::PrintFunctionItem("NN", "nn", nnDist.first);
// #endif
    }

    // free(sizes);
}


// sort the edges by weight when nn is found
// uf save the weight when connected
template<int dim>
inline void naive_link(LDS::EDGE *edges, UnionFind::ParUF<intT> *uf, TreeNNFinder<dim> *finder, intT C){
    //TODO: make parallel
    // find all component ids, find the max
    // sampleSort(edges, C, LDS::edgeComparator());
#ifdef DEBUG
    UTIL::PrintVec2(edges, C);
#endif
    parallel_for(intT i = 1; i < C; ++i){
      // intT u = uf->find(edges[i].first);
      // intT v = uf->find(edges[i].second);
      // if(u==v) continue;
      uf->link(edges[i].first, edges[i].second, edges[i].getW() );
      // intT newc = uf->find(u);
      // finder->merge(u,v,newc);
    }
}

template<int dim>
inline void naive_complete_linkage(TreeNNFinder<dim> *finder, timer t1){
  UTIL::PrintCaption("Naive TREE");
  UnionFind::ParUF<intT> *uf = finder->uf;
  intT n = finder->n;
  LDS::EDGE *edges = (LDS::EDGE *)malloc(sizeof(LDS::EDGE) * n);
  parallel_for(intT i=0; i<n; ++i){
    edges[i] = LDS::EDGE(-1,-1,numeric_limits<double>::max());
  }
#ifdef VERBOSE
	UTIL::PrintSubtimer("initialize", t1.next());
#endif

  int round = 0;
  bool print = false;
  while(finder->C > 1){
    intT C = finder->C;
    round ++;
#ifdef VERBOSE
    if(round < 5 || round % 10 == 1){ //
    print = true;
    UTIL::PrintBreak();
    UTIL::PrintFunctionItem("CLINK", "Comp Num", finder->C);}else{print = false;}
#endif
  if(round == 1){
    finder->findAllNN(edges);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("all-pts-nn", t1.next());
#endif
  }else{
    naive_find_nn_range(edges, finder);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("find-nn", t1.next());
#endif
  }
    naive_link(edges, uf, finder, C);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("link", t1.next());
#endif
    
    finder->updateActiveClustersMerge(round);
#ifdef VERBOSE
	if(print) UTIL::PrintSubtimer("update-active-clusters", t1.next());
#endif
  }
  free(edges);

}
#endif