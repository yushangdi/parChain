#ifndef TEST_H
#define TEST_H

#include "sequence.h"
#include "geometry.h"
#include "sampleSort.h"
#include "shared.h"
#include "unionfind.h"
#include "gettime.h"
#include "kdTree2.h"
#include "neighbor.h"

template<int dim>
inline void kdTreeTestFar(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf){
  timer t1; t1.start();
  typedef iPoint<dim> pointT;
  typedef FINDNN::CLinkNodeInfo nodeInfo;
  typedef kdNode<dim, pointT, nodeInfo> nodeT;
  typedef FINDNN::NNcomplete<dim, nodeT> F;

  F *f = new F();

  pointT *PP = newA(pointT, n);
  parallel_for(intT i=0; i<n; ++i) {
    PP[i] = pointT(P[i], i);}

  // farthest between kdt1 and kdt2
  intT n1 = n/2;
  intT n2 = n-n/2;
  kdTree<dim, pointT, nodeInfo > *kdtree1 = new kdTree<dim, pointT, nodeInfo>(PP, n1);
  kdTree<dim, pointT, nodeInfo > *kdtree2 = new kdTree<dim, pointT, nodeInfo>(PP+n1, n2);
  FINDNN::dualtree_serial<dim, nodeT, F>(kdtree1->root, kdtree2->root, f);

  // farthest within kdt
  // kdTree<dim, pointT, nodeInfo > *kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);
  // FINDNN::dualtree_serial<dim, nodeT, F>(kdtree->root, kdtree->root, f, false);

  LDS::EDGE *e = f->e;
  UTIL::PrintFunctionItem("tree", "x", (*e).first);
  UTIL::PrintFunctionItem("tree", "y", (*e).second);
  UTIL::PrintFunctionItem("tree", "d", (*e).getW());

}
// 2D_GaussianDisc_100K
// [tree] x = 4603
// [tree] y = 59993
// [tree] d = 1741.44
// PBBS-time: 0.0198

template<int dim>
inline void kdTreeTestNear(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf){
  timer t1; t1.start();
  typedef iPoint<dim> pointT;
  typedef FINDNN::CLinkNodeInfo nodeInfo;
  typedef kdNode<dim, pointT, nodeInfo> nodeT;
  typedef FINDNN::AllPtsNN<dim, nodeT> F;
  
  pointT *PP = newA(pointT, n);
  parallel_for(intT i=0; i<n; ++i) {
    PP[i] = pointT(P[i], i);}

  kdTree<dim, pointT, nodeInfo > *kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);

  LDS::EDGE *edges = (LDS::EDGE *)malloc(sizeof(LDS::EDGE) * n);
  parallel_for(intT i=0; i<n; ++i){
    edges[i] = LDS::EDGE(-1,-1,numeric_limits<double>::max());
  }
  F *f = new F(edges);
  FINDNN::dualtree_serial<dim, nodeT, F>(kdtree->root, kdtree->root, f, false);
  for(intT i=0; i<n; ++i) {
   edges[i].print();
  }

}

// kdTreeTestMarkClusters(P,n,uf);
template<int dim>
inline void kdTreeTestMarkClusters(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf){
  typedef iPoint<dim> pointT;
  typedef FINDNN::CLinkNodeInfo nodeInfo;
  typedef kdNode<dim, pointT, nodeInfo> nodeT;

  pointT *PP = newA(pointT, n);
  parallel_for(intT i=0; i<n; ++i) {
    PP[i] = pointT(P[i], i);}

  uf->link(0,1);
  uf->link(5,6);
  uf->link(3,4);

  kdTree<dim, pointT, nodeInfo > *kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);
  FINDNN::markTree<dim, nodeT>(kdtree->root, uf); // mark with point's id
  FINDNN::printTreeCId<dim, nodeT>(kdtree->root);
  FINDNN::markTree<dim, nodeT>(kdtree->root, uf, 4);//mark all with 4
  FINDNN::printTreeCId<dim, nodeT>(kdtree->root);

}

template<int dim>
inline void kdTreeTestNear2Small(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf){
  typedef iPoint<dim> pointT;
  typedef FINDNN::CLinkNodeInfo nodeInfo;
  typedef kdNode<dim, pointT, nodeInfo> nodeT;
  typedef FINDNN::NNsingle<dim, nodeT> F;

  

  pointT *PP = newA(pointT, n);
  parallel_for(intT i=0; i<n; ++i) {
    PP[i] = pointT(P[i], i);}

  pointT *PP2 = newA(pointT, 2);
  PP2[0] = pointT(P[5], 5);
  PP2[1] = pointT(P[6], 6);


  uf->link(0,1);
  uf->link(5,6);
  uf->link(3,4);

  kdTree<dim, pointT, nodeInfo > *kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);
  kdTree<dim, pointT, nodeInfo > *kdtree2 = new kdTree<dim, pointT, nodeInfo>(PP2, 2);

  F *f = new F(uf, uf->find(6));
  FINDNN::markTree<dim, nodeT>(kdtree->root, uf); // mark with point's id
  FINDNN::dualtree_serial<dim, nodeT, F>(kdtree2->root, kdtree->root, f, false);

  LDS::EDGE *e = f->e;
  UTIL::PrintFunctionItem("tree", "x", (*e).first);
  UTIL::PrintFunctionItem("tree", "y", (*e).second);
  UTIL::PrintFunctionItem("tree", "d", (*e).getW());

}

// kdTreeTestNear2(P, n, uf);
template<int dim>
inline void kdTreeTestNear2(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf){
  typedef iPoint<dim> pointT;
  typedef FINDNN::CLinkNodeInfo nodeInfo;
  typedef kdNode<dim, pointT, nodeInfo> nodeT;
  typedef FINDNN::NNsingle<dim, nodeT> F;

  

  pointT *PP = newA(pointT, n);
  parallel_for(intT i=0; i<n; ++i) {
    PP[i] = pointT(P[i], i);}

  pointT *PP2 = newA(pointT,100);
  parallel_for(intT i=0; i<100; ++i) {
    PP2[i] = pointT(P[i], i);
    uf->link(0,i);
  
  }

  kdTree<dim, pointT, nodeInfo > *kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);
  kdTree<dim, pointT, nodeInfo > *kdtree2 = new kdTree<dim, pointT, nodeInfo>(PP2, 100);

  F *f = new F(uf, uf->find(0));
  FINDNN::markTree<dim, nodeT>(kdtree->root, uf); // mark with point's id
  FINDNN::dualtree_serial<dim, nodeT, F>(kdtree2->root, kdtree->root, f, false);

  LDS::EDGE *e = f->e;
  UTIL::PrintFunctionItem("tree", "x", (*e).first);
  UTIL::PrintFunctionItem("tree", "y", (*e).second);
  UTIL::PrintFunctionItem("tree", "d", (*e).getW());

  // pair<pair<intT, intT>, double> result = FINDNN::bruteForceNearest<pointT>(PP2, PP+100, 100, n-100);
  // UTIL::PrintFunctionItem("tree", "x", result.first.first);
  // UTIL::PrintFunctionItem("tree", "y", result.first.second);
  // UTIL::PrintFunctionItem("tree", "d", result.second);
  // [tree] x = 97
  // [tree] y = 9106
  // [tree] d = 0.391903
}

template<int dim>
inline void nnCandidateTestSmall(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf){
  typedef iPoint<dim> pointT;
  typedef FINDNN::CLinkNodeInfo nodeInfo;
  typedef kdNode<dim, pointT, nodeInfo> nodeT;

  pointT *PP = newA(pointT, n);
  parallel_for(intT i=0; i<n; ++i) {
    PP[i] = pointT(P[i], i);}

  uf->link(0,1);
  uf->link(5,6);
  uf->link(3,4);

  kdTree<dim, pointT, nodeInfo > *kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);
  FINDNN::markTree<dim, nodeT>(kdtree->root, uf); // mark with point's id

  intT *arr = newA(intT, n);
  typedef FINDNN::RangeQueryF<pointT, nodeT> Fr;
  Fr *fr = new Fr(uf, arr, n);
  FINDNN::NNcandidate<dim, pointT, Fr, nodeT>(kdtree->root, PP+2, 3, fr);

  UTIL::PrintVec(fr->sizes, n);


}

template<int dim>
inline void nnCandidateTest(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf){
  typedef iPoint<dim> pointT;
  typedef FINDNN::CLinkNodeInfo nodeInfo;
  typedef kdNode<dim, pointT, nodeInfo> nodeT;

  pointT *PP = newA(pointT, n);
  parallel_for(intT i=0; i<n; ++i) {
    PP[i] = pointT(P[i], i);}

  uf->link(2253,4420);
  uf->link(5401,7412);
  uf->link(5401,8652);
  uf->link(5401,1);

  //r==50:
  // 2253
  // 4420
  // 5401
  // 6679
  // 7412
  // 8652

  double r = 100;

  kdTree<dim, pointT, nodeInfo > *kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);
  FINDNN::markTree<dim, nodeT>(kdtree->root, uf); // mark with point's id

  intT *arr = newA(intT, n);
  typedef FINDNN::RangeQueryF<pointT, nodeT> Fr;
  Fr *fr = new Fr(uf, arr, n);
  FINDNN::NNcandidate<dim, pointT, Fr, nodeT>(kdtree->root, PP, r, fr);

  intT* true_sizes = FINDNN::bruteForceNNCandidateBox(P[0], P, n, r, uf);

  parallel_for(intT i=0; i<n; ++i){
    if(true_sizes[i] != arr[i]) cout << i << " " << true_sizes[i]<<" " << arr[i] <<endl;
  }

  free(arr);
  free(true_sizes);

}


void testSerializeUF(UnionFind::ParUF<intT> *uf){
  uf->link(1,2);
  uf->serialize();
  UnionFind::ParUF<intT> *uf2 = new UnionFind::ParUF<intT>(uf->m_n, true);
  uf2->deserialize();
}


// testRangeBB(finder);
template<int dim>
void testRangeBB(TreeNNFinder<dim> *finder){
    typedef iPoint<dim> pointT;
    typedef FINDNN::CLinkNodeInfo nodeInfo;
    typedef kdNode<dim, iPoint<dim>, nodeInfo> nodeT;

    typedef FINDNN::MarkClusterId<dim, nodeT> M;
    M *cidMarker = new M(finder->uf);
    FINDNN::singletree<nodeT, M, intT>(finder->kdtree->root, cidMarker, -1);

    finder->uf->link(0,1);
    finder->uf->link(3,4);
    finder->uf->link(5,6);

    pointT *items = newA(pointT, 3);
    LDS::EDGE *edges = newA(LDS::EDGE, 3);
    intT *CMap = newA(intT, 7);

    intT cid = finder->uf->find(0);
    items[0] = (finder->PP[cid]);
    edges[0] = LDS::EDGE(cid, 2, 1.45);
    CMap[cid] = 0;

    cid = finder->uf->find(3);
    items[1] = (finder->PP[cid]);
    edges[1] = LDS::EDGE(cid, 2, 4.5);
    CMap[cid] = 1;

    cid = finder->uf->find(5);
    items[2] = (finder->PP[cid]);
    edges[2] = LDS::EDGE(cid, 2, 5.5);
    CMap[cid] = 2;


    kdTree<dim, pointT, nodeInfo>* queryTree = new kdTree<dim, pointT, nodeInfo>(items, 3);
    typedef FINDNN::RangeBoundingBox<dim, nodeT> M2;
    M2 *Mbb = new M2(edges, CMap);
    FINDNN::singletree<nodeT, M2, intT>(queryTree->root, Mbb, -1);

    FINDNN::printTreeCId<dim, nodeT>(queryTree->root);
}


  // allNNDebug2(P, n, uf, true);
  // exit(1);
template<int dim>
void allNNDebug2(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, bool par = true) {
    cout << "num workers " << getWorkers() << endl;
	typedef FINDNN::CLinkNodeInfo nodeInfo;
  typedef iPoint<dim> pointT;
  typedef kdNode<dim, pointT, nodeInfo> nodeT;

  //copy the points into PP array
  pointT *PP = newA(pointT, n);
  parallel_for(intT i=0; i<n; ++i) {
      PP[i] = pointT(P[i], i);
  }

  LDS::EDGE *edges = newA(LDS::EDGE, n); // the nn of point with id i is stored at edges[i]
  parallel_for(intT i=0; i<n; ++i){
    edges[i] = LDS::EDGE(-1,-1,numeric_limits<double>::max());
  }
  kdTree<dim, pointT, nodeInfo> *kdtree = new kdTree<dim, pointT, nodeInfo>(PP, n);

intT *corr = newA(intT, n);
  typedef FINDNN::AllPtsNN<dim, nodeT> F;
  F *f = new F(edges); // store results in edges array
  if(par){
std::ifstream file{"edges_correct_5M.txt"};
 intT u, v;
 size_t line_ct = 0;
    while (file >> u >> v) {
    corr[u] = v;
    if(line_ct != u) cout << "!!! " << endl;
    line_ct ++;
    }
cout << "loaded edges" << endl;

    FINDNNP::dualtree<dim, nodeT>(kdtree->root, kdtree->root, f);
cout << "found nn" << endl;
    
LDS::edgeComparator2 EC2 = LDS::edgeComparator2();
    line_ct = 0;
    while (line_ct < n) {
	    u = line_ct;v = corr[u];
      if(u != edges[line_ct].first || v != edges[line_ct].second){
        cout << "=====" << endl;
        cout << u << " " << v << " " << P[u].pointDist(P[v]) <<  endl;
        cout << edges[line_ct].first << " " <<  edges[line_ct].second<< " " << edges[line_ct].getW() << " "  << P[u].pointDist(P[edges[line_ct].second]) << endl;
      }
      line_ct ++;
    }

  }else{
    FINDNN::dualtree_serial<dim, nodeT>(kdtree->root, kdtree->root, f);
    std::ofstream out{"edges_correct_5M.txt"};
    for (intT i = 0; i < n; ++i) {
      out << edges[i].first << " " << edges[i].second << std::endl;
    } 
    out.close();
    std::cout << "# Wrote file." << std::endl;
    }

  free(edges);

}

// using avg.pbbs
template<int dim>
void testAverageDistance(point<dim>* P, intT n){
  typedef iPoint<dim> pointT;
  typedef FINDNN::AveLinkNodeInfo<dim, pointT *> nodeT;

  pointT *PP = newA(pointT, n);
  for(intT i=0; i<n; ++i) {PP[i] = pointT(P[i], i);}
  LDS::node<pointT *> *PPLList = newA(LDS::node<pointT *>, n);
  for(intT i=0; i<n; ++i) {PPLList[i] = LDS::node<pointT *>(PP+i, nullptr);}
  // PPLList[0].next = PPLList + 4;
  // PPLList[4].next = PPLList + 5;
  // PPLList[1].next = PPLList + 3;
  // PPLList[2].next = PPLList + 1;


  nodeT *nodes = newA(nodeT, 2*n);
  parallel_for(intT i=0; i<n; ++i) { nodes[i] = nodeT(i, PPLList+i);}
  nodes[6] = nodeT(3, 1, 6, nodes, nodes + 4);
  nodes[7] = nodeT(4, 2, 7, nodes+6, nodes + 5);
  nodes[8] = nodeT(8, 1, 8, nodes+1, nodes + 3);
  nodes[9] = nodeT(8, 2, 9, nodes+2, nodes + 8);
  // parallel_for(intT i=0; i<n; ++i) { FINDNN::ConstructAveLinkNodeInfo<pointT *>(i, nodes+i, PPLList+i);}
  // FINDNN::ConstructAveLinkNodeInfo<pointT *>(3, 1, 6, nodes+6, nodes, nodes + 4);
  // FINDNN::ConstructAveLinkNodeInfo<pointT *>(4, 2, 7, nodes+7, nodes+6, nodes + 5);
  // FINDNN::ConstructAveLinkNodeInfo<pointT *>(8, 1, 8, nodes+8, nodes+1, nodes + 3);
  // FINDNN::ConstructAveLinkNodeInfo<pointT *>(8, 2, 9, nodes+9, nodes+2, nodes + 8);
  // for(intT i=0; i<10; ++i) {cout << (nodes+i)->items->elt->x(1) << endl;}

  nodeT *inode = nodes + 7;
  nodeT *jnode = nodes + 9;
  typedef FINDNN::NNaverage<dim, nodeT, pointT *> F;
  F *f = new F();
  double result = FINDNN::dualtree_serial2<nodeT, F, double>(inode, jnode, f);
  cout << result << endl;
}

//using avg2.pbbs
template<int dim>
void testAveRangeQuerry(point<dim>* P, intT n, TreeNNFinder<dim> *finder){
  typedef iPoint<dim> pointT;
  typedef FINDNN::AveLinkNodeInfo<dim, pointT *> nodeT;
  typedef kdNode<dim, pointT, FINDNN::CLinkNodeInfo> kdnodeT;

  pointT *PP = finder->PP;
  LDS::node<pointT *> *PPLList = finder->PPLList;

  nodeT *nodes = finder->nodes;
  nodes[6] = nodeT(1, 1, 6, nodes, nodes + 1);
  nodes[7] = nodeT(3, 1, 7, nodes+2, nodes + 3);
  nodes[8] = nodeT(5, 1, 8, nodes+4, nodes + 5);

  finder->uf->link(0,1);
  finder->uf->link(2,3);
  finder->uf->link(4,5);

  finder->rootIdx[1] = 6;
  finder->rootIdx[3] = 7;
  finder->rootIdx[5] = 8;

  finder->cacheTbs[6] = new LDS::distCacheT(3,LDS::hashClusterAve());
  finder->cacheTbs[7] = new LDS::distCacheT(3,LDS::hashClusterAve());
  finder->cacheTbs[8] = new LDS::distCacheT(3,LDS::hashClusterAve());

  typedef FINDNN::RangeQueryHashAveF<dim, pointT, kdnodeT, nodeT> Fr;
  Fr *fr = new Fr(finder->uf, 1, nodes, finder->rootIdx, finder->cacheTbs, finder->edges); 
  FINDNN::NNcandidate<dim, Fr, kdnodeT, nodeT>(finder->kdtree->root, nodes + 6, 3, fr);
  fr = new Fr(finder->uf, 3, nodes, finder->rootIdx, finder->cacheTbs, finder->edges); 
  FINDNN::NNcandidate<dim, Fr, kdnodeT, nodeT>(finder->kdtree->root, nodes + 7, 4, fr);
  fr = new Fr(finder->uf, 5, nodes, finder->rootIdx, finder->cacheTbs, finder->edges); 
  FINDNN::NNcandidate<dim, Fr, kdnodeT, nodeT>(finder->kdtree->root, nodes + 8, 2.5, fr);

  UTIL::PrintVec2<LDS::EDGE>(finder->edges, n);
  cout << "==" << endl;
  UTIL::PrintTable<LDS::distCacheT>(finder->cacheTbs[6]);
  cout << "==" << endl;
  UTIL::PrintTable<LDS::distCacheT>(finder->cacheTbs[7]);
  cout << "==" << endl;
  UTIL::PrintTable<LDS::distCacheT>(finder->cacheTbs[8]);

}

template< class TF>
inline void testTableFull(TF *finder){
  intT count = 0;
  intT C = finder->C;
  if(finder->C > 50) return;
  for(intT i=0; i<C; ++i){
    intT cid1 = finder->activeClusters[i];
    for(intT j = i+1; j < C; ++j){
      intT cid2 = finder->activeClusters[j];
      double d = finder->find(cid1, cid2);
      if(d == UNFOUND_TOKEN || d == CHECK_TOKEN){
        if(count < 10) cout << cid1 << " " << cid2 << endl;
        count++;
      }
    }
  }
}


#endif