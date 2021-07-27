#include "sequence.h"
#include "geometry.h"
#include "sampleSort.h"
#include "parseCommandLine.h"

#define VERBOSE
// #define DEBUG
#define BENCHCACHE
#define PERF_RANGE
// #define ELTPERCACHELINE 128/sizeof(intT)

#define ELTPERCACHELINE sizeof(intT)/sizeof(intT)

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
// #include "dynamicKdTree.h"
#include "neighbor.h"
// #include "distanceMatrix.h"
#include "method_chain_tree.h"
#include "method_chain_matrix_range.h"
// #include "method_boruvka.h"
#include "dendrogram.h"

#include <limits>
#include "method_chain_matrix.h"

using namespace std;


// *************************************************************
//    DRIVER
// *************************************************************



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
        ChainMatrix::completeLinkageMatrix<dim, pointT>(P, n, uf, eps);
    }else if(method ==  "ward"){
        ChainMatrix::wardLinkageMatrix<dim, pointT>(P, n, uf, eps);
    }else if(method ==  "avg"){
        ChainMatrix::avgLinkageMatrix<dim, pointT>(P, n, uf, eps);
    }else if(method ==  "avgsq"){
        ChainMatrix::avgsqLinkageMatrix<dim, pointT>(P, n, uf, eps);
    }else{
        cout << "invalid method" << endl;
        exit(1);
    }
    cout << std::fixed;
    cout << std::setprecision(10);
    UTIL::PrintFunctionItem("CLINK", "COST", uf->cost() );
    return uf;
  }

  bool use_matrix_range = params.getOption("-matrixrange");
  cout << "use matrix range = " <<  use_matrix_range << endl;

  intT naive_thresh = params.getOptionIntValue("-naivethresh", 5);
  if(naive_thresh == 1){ //0: remove range query optimization
    naive_thresh = n;
  }
  cout << "naive_thresh = " <<  naive_thresh << endl;

 if(use_matrix_range){
    if(method ==  "complete"){
        ChainMatrixRange::completeLinkage<dim, pointT>(P, n, uf, eps, naive_thresh);
    }else if(method ==  "ward"){
        ChainMatrixRange::wardLinkage<dim, pointT>(P, n, uf, eps, naive_thresh);
    }else if(method ==  "avg"){
        ChainMatrixRange::avgLinkage<dim, pointT>(P, n, uf, eps, naive_thresh);
    }else if(method ==  "avgsq"){
        ChainMatrixRange::avgsqLinkage<dim, pointT>(P, n, uf, eps, naive_thresh);
    }else{
        cout << "invalid method" << endl;
        exit(1);
    }
    cout << std::fixed;
    cout << std::setprecision(10);
    UTIL::PrintFunctionItem("CLINK", "COST", uf->cost() );
    return uf;
}


  intT cache_size = params.getOptionIntValue("-cachesize", 32);
  if(cache_size == 1) cache_size =0;
  cout << "cache_size/2 = " <<  cache_size << endl;

  if(method ==  "complete"){
      ChainTree::completeLinkage<dim, pointT>(P, n, uf, eps, naive_thresh, cache_size);
  }else if(method ==  "ward"){
      ChainTree::wardLinkage<dim, pointT>(P, n, uf, eps, naive_thresh, cache_size);
  }else if(method ==  "avg"){
      ChainTree::avgLinkage<dim, pointT>(P, n, uf, eps, naive_thresh, cache_size);
  }else if(method ==  "avgsq"){
      ChainTree::avgsqLinkage<dim, pointT>(P, n, uf, eps, naive_thresh, cache_size);
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






