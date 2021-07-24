
#include <stdlib.h>  
#include <stdio.h>  
#include <limits>
#include <tuple>
#include <queue>
#include <vector>
#include <map>
#include "util.h"
#include "sampleSort.h"
#include "unionFind.h"
#include "gettime.h"
using namespace std;

// https://github.com/scikit-learn/scikit-learn/blob/2beed55847ee70d363bdbfe14ee4401438fba057/sklearn/cluster/_agglomerative.py#L346
// https://github.com/scikit-learn/scikit-learn/blob/2beed55847ee70d363bdbfe14ee4401438fba057/sklearn/cluster/_hierarchical_fast.pyx

// to build a complete tree, sklearn call scipy for ward's linkage

template<class F>
map<int, double> *merge(map<int, double> *a, map<int, double> *b, int* mask, int n_a, int n_b, F& new_dist ){
  map<int, double> *out_obj = new map<int, double>();
  for(auto iter = a->begin(); iter != a->end(); ++iter){
    int key = iter->first;
    if(mask[key]>0){
      out_obj->insert_or_assign({key, iter->second});
    }
  }

  for(auto iter = b->begin(); iter != b->end(); ++iter){
    int key = iter->first;
    double value = iter->second;
    if(mask[key]>0){
      auto out_it = out_obj->find(key);
      if(out_it == b->end()){
        out_obj->insert({key, value});
      }else{
        out_obj->insert_or_assign({key, new_dist(out_it->second, value, 0, n_a, n_b, 0)});
      }
      
    }
  }

  return out_obj;
}


template<class F>
void sklearn(double * D, int n, F new_dist){
  typedef tuple<double, int, int> entry;
  int n_samples = n;
  int n_ndoes = 2*n_samples-1;
  map<int, double> **A = newA((map<int, double> *), n_nodes);
  priority_queue<entry> inertia;
  for (intT i=0; i<n; i++) {
    A[i] = new map<int, double>();
    for (intT j=i+1; j<n; j++) {
      double d = D[condensed_index(n,i,j)];
      A[i]->insert({j,d});
      inertia.push(entry(d, i, j));
    }
  }

  vector<pair<int, int>> children;
  int* parent = newA(int, n_ndoes);
  int* used_node = newA(int, n_ndoes);
  parallel_for(int i=0; i<n_ndoes; ++i){
    used_node[i] = 1;
    parent[i] = i;
  }

  int i,j = 0;
  int n_i,n_j = 1;
  entry edge;
  double cost = 0;
  for(int k=n_samples;k<n_nodes; ++k){
    //identify the merge
    while(true){
      edge = inertia.pop();
      if(used_node[get<1>(edge)] > 0 && used_node[get<2>(edge)] > 0) break
    }
    i = get<1>(edge);
    j = get<2>(edge);
    cost += get<0>(edge);

    parent[i]=k;
    parent[j]=k;

    children.emplace_back(make_pair(i,j));

    //Keep track of the number of elements per cluster
    n_i = used_node[i];
    n_j = used_node[j];
    used_node[k] = n_i + n_j;
    used_node[i] = used_node[j] = 0;

    // update the structure matrix A and the inertia matrix
    // a clever 'min', or 'max' operation between A[i] and A[j]
    map<int, double> *coord_col = merge(A[i], A[j], used_node, n_i, n_j, new_dist);
    for(auto iter = coord_col->begin(); iter != coord_col->end(); ++iter){
      int col = iter->first;
      double d = iter->second;
      A[col]->insert_or_assign({k,d});
      inertia.push(entry(d,k,col));
    }
    A[k] = coord_col;
    delete A[i];
    delete A[j];
  }
  // int n_leaves = n_samples;
  free(parent);
  free(used_node);
  free(A);

}