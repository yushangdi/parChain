#include "sequence.h"
#include "geometry.h"
#include "sampleSort.h"
#include "parseCommandLine.h"

#define VERBOSE
// #define DEBUG

#ifdef USEJEMALLOC
#include<jemalloc/jemalloc.h>
#define jeNewA(__E,__n) (__E*) je_custom_prefix_malloc((__n)*sizeof(__E))
#define jeFree(__E) je_custom_prefix_free(__E)
#endif


#define PRINT_FREQ 1
#define LINKAGE_DOPRINT(round) (round < 5 || round % PRINT_FREQ == 0)
#define INTT_MAX numeric_limits<intT>::max()
#define TIMING2

#include "shared.h"
#include "unionfind.h"
#include "gettime.h"
#include "kdTree2.h"
#include "dynamicKdTree.h"
#include "neighbor.h"
// #include "distanceMatrix.h"
#include "method_chain_tree.h"
// #include "method_boruvka.h"
#include "dendrogram.h"

#include <limits>
#include "method_chain_matrix.h"

using namespace std;


// *************************************************************
//    DRIVER
// *************************************************************

template<int dim, class pointT>
inline void completeLinkage(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps, intT naive_thresh, intT cache_size){
  timer t1;
  t1.start();
  cout << "complete Linkage of " << n << ", dim " << dim << " points" << endl; 
  using nodeInfo = FINDNN::CLinkNodeInfo; // for kdtree
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>; // point array
  using distT = distComplete3<dim, pointT, nodeT, nodeInfo>; //kdtree
  using boxT = FINDNN::queryBallSimple<dim, nodeT>;
  bool no_cache = true;
  using Fr = FINDNN::RangeQueryCountF1<dim, pointT, nodeInfo, distT, boxT>;
  using M = FINDNN::MarkClusterId<dim, Fr>;

  using FinderT = TreeNNFinder<dim, distT, Fr, M>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, no_cache, eps, naive_thresh, cache_size); // if true, do not keep cache
  ChainTree::chain_linkage<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());
}

template<int dim, class pointT>
inline void wardLinkage(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps, intT naive_thresh, intT cache_size){
  timer t1;
  t1.start();
  cout << "ward Linkage of " << n << ", dim " << dim << " points" << endl;
  using nodeInfo = FINDNN::WLinkNodeInfo; // for kdtree
  using nodeT = FINDNN::LinkNodeInfo3<dim, pointT *>;
  using MCenter = FINDNN::MarkKdTreeCenters<dim, pointT, nodeInfo>;
  using distT = distWard1<dim, pointT, nodeT, nodeInfo, MCenter>;
  using boxT = FINDNN::queryBallWard<dim, nodeT>;
  bool no_cache = true;
  using Fr = FINDNN::RangeQueryCenterF<dim, pointT, nodeInfo, distT, boxT>;
  using M = FINDNN::DummyMarker<Fr>;

  using FinderT = TreeNNFinder<dim, distT, Fr, M>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, no_cache, eps, naive_thresh, cache_size); // if true, do not keep cache
  ChainTree::chain_linkage<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());
}

template<int dim, class pointT>
inline void avgLinkage(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps, intT naive_thresh, intT cache_size){
  timer t1;
  t1.start();
  cout << "average Linkage of " << n << ", dim " << dim << " points" << endl; 
  using nodeInfo = FINDNN::CLinkNodeInfo; // for kdtree
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;
  using distT = distAverage5<dim, pointT, nodeT>; // consecutive point array
  using boxT = FINDNN::queryBallSimple<dim, nodeT>;
  bool no_cache = cache_size == 0;
  using Fr = FINDNN::RangeQueryCenterF<dim, pointT, nodeInfo, distT, boxT>;
  using M = FINDNN::DummyMarker<Fr>;

  using FinderT = TreeNNFinder<dim, distT, Fr, M>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, no_cache, eps, naive_thresh, cache_size); // if true, do not keep cache
  ChainTree::chain_linkage<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());
}

template<int dim, class pointT>
inline void avgsqLinkage(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps, intT naive_thresh, intT cache_size){
  timer t1;
  t1.start();
  cout << "average Linkage of " << n << ", dim " << dim << " points, sqeuclidean" << endl; 
  using nodeInfo = FINDNN::CLinkNodeInfo; // for kdtree
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;
  using distT = distAverage4<dim, pointT, nodeT, nodeInfo>; // consecutive point array
  using boxT = FINDNN::queryBallSimple<dim, nodeT>; //still use simple, distT will  postprocess
  bool no_cache = true;
  using Fr = FINDNN::RangeQueryCenterF<dim, pointT, nodeInfo, distT, boxT>;
  using M = FINDNN::DummyMarker<Fr>;

  using FinderT = TreeNNFinder<dim, distT, Fr, M>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, no_cache, eps, naive_thresh, cache_size); // if true, do not keep cache
  ChainTree::chain_linkage<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());
}

template<int dim, class pointT>
inline void completeLinkageMatrix(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps){
  timer t1;
  t1.start();
    // DISTANCE MATRIX
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;

  // cout << "Matrix Distance Average Linkage of " << n << ", dim " << dim << " points" << endl; 
  // using distT = distAverage1<dim, pointT, nodeT>;

  cout << "Matrix Distance Complete Linkage of " << n << ", dim " << dim << " points" << endl; 
  using distT = distComplete3<dim, pointT, nodeT, FINDNN::dummyNodeInfo>;

  // cout << "Matrix Distance Ward Linkage of " << n << ", dim " << dim << " points" << endl;
  // using distT = distWard1<dim, pointT, nodeT, FINDNN::dummyNodeInfo, FINDNN::DummyMarkerCenters>;

  using FinderT = MatrixNNFinder<dim, distT>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, false, eps); 
  ChainMatrix::chain_linkage_matrix<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());

}

template<int dim, class pointT>
inline void wardLinkageMatrix(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps){
  timer t1;
  t1.start();
    // DISTANCE MATRIX
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;

  // cout << "Matrix Distance Average Linkage of " << n << ", dim " << dim << " points" << endl; 
  // using distT = distAverage1<dim, pointT, nodeT>;

  cout << "Matrix Distance Ward Linkage of " << n << ", dim " << dim << " points" << endl;
  using distT = distWard1<dim, pointT, nodeT, FINDNN::dummyNodeInfo, FINDNN::DummyMarkerCenters>;

  using FinderT = MatrixNNFinder<dim, distT>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, false, eps); 
  ChainMatrix::chain_linkage_matrix<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());

}

template<int dim, class pointT>
inline void avgLinkageMatrix(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps){
  timer t1;
  t1.start();
    // DISTANCE MATRIX
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;

  cout << "Matrix Distance Average Square Linkage of " << n << ", dim " << dim << " points" << endl; 
  // using distT = distAverage1<dim, pointT, nodeT>;
  using distT = distAverage5<dim, pointT, nodeT>;

  using FinderT = MatrixNNFinder<dim, distT>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, false, eps); 
  ChainMatrix::chain_linkage_matrix<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());

}

template<int dim, class pointT>
inline void avgsqLinkageMatrix(point<dim>* P, intT n, UnionFind::ParUF<intT> *uf, double eps){
  timer t1;
  t1.start();
    // DISTANCE MATRIX
  using nodeT = FINDNN::LinkNodeInfo1<dim, pointT *>;

  cout << "Matrix Distance Average Square Linkage of " << n << ", dim " << dim << " points" << endl; 
  using nodeInfo = FINDNN::CLinkNodeInfo; // for kdtree
  using distT = distAverage4<dim, pointT, nodeT, nodeInfo>; // consecutive point array

  using FinderT = MatrixNNFinder<dim, distT>;
  distT::printName();
  FinderT *finder = new FinderT(n, P, uf, false, eps); 
  ChainMatrix::chain_linkage_matrix<dim, FinderT>(finder, t1);
  timer t2;t2.start();
  dendrogram::dendroLine* dendro = dendrogram::formatDendrogram<nodeT>(finder->nodes, n, eps);
  UTIL::PrintSubtimer("format-dendro", t2.next());
  free(dendro);
  delete finder;
  UTIL::PrintSubtimer("CLINK", t1.next());

}

template<int dim>
UnionFind::ParUF<intT> *linkage(point<dim>* P, intT n, commandLine params, UnionFind::ParUF<intT> *uf) {
  cout << "+++++++++++++++" << endl;
  if (n < 2) abort();
  typedef iPoint<dim> pointT;

  double eps = params.getOptionDoubleValue("-eps", 1e-20);
  cout << "eps = " << std::setprecision(25) <<  eps << endl;

  bool use_matrix = params.getOption("-matrix");
  cout << "use matrix = " <<  use_matrix << endl;

  string method = params.getOptionValue("-method", "invalid");

  if(use_matrix){
    if(method ==  "complete"){
        completeLinkageMatrix<dim, pointT>(P, n, uf, eps);
    }else if(method ==  "ward"){
        wardLinkageMatrix<dim, pointT>(P, n, uf, eps);
    }else if(method ==  "avg"){
        avgLinkageMatrix<dim, pointT>(P, n, uf, eps);
    }else if(method ==  "avgsq"){
        avgsqLinkageMatrix<dim, pointT>(P, n, uf, eps);
    }else{
        cout << "invalid method" << endl;
        exit(1);
    }
    cout << std::fixed;
    cout << std::setprecision(10);
    UTIL::PrintFunctionItem("CLINK", "COST", uf->cost() );
    return uf;
  }

  intT naive_thresh = params.getOptionIntValue("-naivethresh", 5);
  if(naive_thresh == 1){ //0: remove range query optimization
    naive_thresh = n;
  }
  cout << "naive_thresh = " <<  naive_thresh << endl;

  intT cache_size = params.getOptionIntValue("-cachesize", 32);
  if(cache_size == 1) cache_size =0;
  cout << "cache_size = " <<  cache_size << endl;

  if(method ==  "complete"){
      completeLinkage<dim, pointT>(P, n, uf, eps, naive_thresh, cache_size);
  }else if(method ==  "ward"){
      wardLinkage<dim, pointT>(P, n, uf, eps, naive_thresh, cache_size);
  }else if(method ==  "avg"){
      avgLinkage<dim, pointT>(P, n, uf, eps, naive_thresh, cache_size);
  }else if(method ==  "avgsq"){
      avgsqLinkage<dim, pointT>(P, n, uf, eps, naive_thresh, cache_size);
  }else{
      cout << "invalid method" << endl;
      exit(1);
  }

  // // uf->serialize("./uf_" + version + ".txt" );
  cout << std::fixed;
  cout << std::setprecision(10);
  UTIL::PrintFunctionItem("CLINK", "COST", uf->cost() );

  // uf->deserialize("./uf_" + version + ".txt" );
  // intT cutoff = 10000;
  // uf->fcluster(cutoff, true, "./cluster_" + version + "_"  + std::to_string(cutoff) + ".txt");
  return uf;
}

template UnionFind::ParUF<int> *linkage<2>(point<2>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<3>(point<3>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<4>(point<4>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<5>(point<5>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<6>(point<6>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<7>(point<7>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<8>(point<8>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<9>(point<9>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<10>(point<10>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<16>(point<16>*, intT, commandLine, UnionFind::ParUF<intT> *);
template UnionFind::ParUF<int> *linkage<128>(point<128>*, intT, commandLine, UnionFind::ParUF<intT> *);






