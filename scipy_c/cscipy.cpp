
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

#define COMP 1
#define WARD 2
#define AVG 3
#define AVGSQ 4



// https://github.com/scipy/scipy/blob/v1.7.0/scipy/cluster/_hierarchy.pyx
// https://github.com/scipy/scipy/blob/701ffcc8a6f04509d115aac5e5681c538b5265a2/scipy/cluster/_hierarchy_distance_update.pxi

#define NPY_INFINITYF numeric_limits<double>::max()

typedef tuple<int, int, double, int> dendroentry;

void label(dendroentry* Z, int n){
    //"""Correctly label clusters in unsorted dendrogram."""
    unionFind uf = unionFind(n);
    int x, y, x_root, y_root = 0;
    for(int i=0; i<n-1; ++i){
        x = get<0>(Z[i]); y = get<1>(Z[i]);
        x_root = uf.find(x); y_root =  uf.find(y);
        if (x_root < y_root){
            get<0>(Z[i]) = x_root; get<1>(Z[i]) = y_root;
        }else{
            get<0>(Z[i]) = y_root;  get<1>(Z[i]) = x_root;
        }
        uf.link(x_root, y_root);
        get<3>(Z[i]) = uf.find(x_root);
    }
}


int condensed_index(int n, int i, int j){
    if (i < j)return n * i - (i * (i + 1) / 2) + (j - i - 1);
    if( i > j)return n * j - (j * (j + 1) / 2) + (i - j - 1);
    cout << "i==j in indexing" << endl;
    return -1;
}
        

// """
// dist will be modified
// """
template<class F>
dendroentry* nn_chain(double * D, int n, F new_dist){
  
  dendroentry *Z_arr = newA(dendroentry, (n-1));
  int *size = new int[n];
  int *cluster_chain = new int[n];
  parallel_for(int i=0; i<n; ++i){
    size[i] = 1;
    cluster_chain[i] = 0;
  }
  int chain_length = 0;
  int x, y, nx, ny, ni = 0;
  double dist, current_min = 0;
  double cost = 0;

  for(int k=0; k<n-1; ++k){

    // get the start of the chain
    if(chain_length == 0){
      chain_length = 1;
      for(int i=0; i<n; ++i){
        if(size[i] > 0){
          cluster_chain[0]=i;
          break;
        }
      }
    }

    // Go through chain of neighbors until two mutual neighbors are found.
    while(true){
      x = cluster_chain[chain_length-1];

      // We want to prefer the previous element in the chain as the
      // minimum, to avoid potentially going in cycles.
      if (chain_length > 1){
        y = cluster_chain[chain_length - 2];
        current_min = D[condensed_index(n, x, y)];
      }else{
        current_min = NPY_INFINITYF;
      }

      for(int i=0; i<n; ++i){
        if(size[i] == 0 || x == i)continue;
        
        dist = D[condensed_index(n, x, i)];
        if(dist < current_min){
          current_min = dist;
          y=i;
        }
      }

      if (chain_length > 1 && y == cluster_chain[chain_length -2])break;

      cluster_chain[chain_length] = y;
      chain_length += 1;

    }//end while true

    //Merge clusters x and y and pop them from stack.
    chain_length -= 2;

    //This is a convention used in fastcluster.
    if(x > y)swap(x,y);

    // get the original numbers of points in clusters x and y
    nx = size[x];
    ny = size[y];

    // Record the new node.
    Z_arr[k] = dendroentry(x,y,current_min, nx+ny);
    size[x] = 0;  // Cluster x will be dropped.
    size[y] = nx + ny;  // Cluster y will be replaced with the new cluster
    // cout << x << " " << y << " " << current_min << endl;
    cost += current_min;

    //Update the distance matrix.
    for(int i=0; i<n; ++i){
      ni = size[i];
      if(ni == 0 || i==y) continue;
      D[condensed_index(n, i, y)] = new_dist(D[condensed_index(n, i, x)],D[condensed_index(n, i, y)],current_min, nx, ny, ni);
    }



  }//end n-1 loop for cluster

  cout << "cost: "  << std::setprecision(10) << cost << endl;
  // Sort Z by cluster distances.
  // order = np.argsort(Z_arr[:, 2], kind='mergesort')
  // Z_arr = Z_arr[order]
  sampleSort(Z_arr, n-1, [&](dendroentry &i, dendroentry &j){return get<2>(i) < get<2>(j);} );

  // Find correct cluster labels inplace.
  label(Z_arr, n);

  free(size);
  free(cluster_chain);
  return Z_arr;

}



struct Complete{
  double operator ()(double d_xi, double d_yi, double d_xy, int nx, int ny, int ni){
    return max(d_xi, d_yi);
  }
};


struct Ward{
  double operator ()(double d_xi, double d_yi, double d_xy, int size_x, int size_y, int size_i){
    double t = 1.0 / (size_x + size_y + size_i);
    return sqrt((size_i + size_x) * t * d_xi * d_xi +
                (size_i + size_y) * t * d_yi * d_yi -
                size_i * t * d_xy * d_xy);
  }
};


struct Average{
  double operator ()(double d_xi, double d_yi, double d_xy, int nx, int ny, int ni){
    return (nx * d_xi + ny * d_yi) / (nx + ny);
  }
};




/////////////////////////////////////// SKLEARN

// https://github.com/scikit-learn/scikit-learn/blob/2beed55847ee70d363bdbfe14ee4401438fba057/sklearn/cluster/_agglomerative.py#L346
// https://github.com/scikit-learn/scikit-learn/blob/2beed55847ee70d363bdbfe14ee4401438fba057/sklearn/cluster/_hierarchical_fast.pyx

// to build a complete tree, sklearn call scipy for ward's linkage

template<class F>
map<int, double> *merge(map<int, double> *a, map<int, double> *b, int* mask, int n_a, int n_b, F& new_dist ){
  map<int, double> *out_obj = new map<int, double>();
  for(auto iter = a->begin(); iter != a->end(); ++iter){
    int key = iter->first;
    if(mask[key]>0){
      out_obj->insert_or_assign(key, iter->second);
    }
  }

  for(auto iter = b->begin(); iter != b->end(); ++iter){
    int key = iter->first;
    double value = iter->second;
    if(mask[key]>0){
      auto out_it = out_obj->find(key);
      if(out_it == b->end()){
        out_obj->insert_or_assign(key, value);
      }else{
        out_obj->insert_or_assign(key, new_dist(out_it->second, value, 0, n_a, n_b, 0));
      }
      
    }
  }

  return out_obj;
}


template<class F>
void sklearn(double * D, int n, F new_dist){
  typedef tuple<double, int, int> entry;
  int n_samples = n;
  int n_nodes = 2*n_samples-1;
  map<int, double> **A = (map<int, double> **) malloc((n_nodes)*sizeof(map<int, double> *));
  priority_queue<entry> inertia;
  for (intT i=0; i<n; i++) {
    A[i] = new map<int, double>();
    for (intT j=i+1; j<n; j++) {
      double d = D[condensed_index(n,i,j)];
      A[i]->insert_or_assign(j,d);
      inertia.push(entry(d, i, j));
    }
  }

  vector<pair<int, int>> children;
  int* parent = newA(int, n_nodes);
  int* used_node = newA(int, n_nodes);
  parallel_for(int i=0; i<n_nodes; ++i){
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
      edge = inertia.top();
      inertia.pop();
      if(used_node[get<1>(edge)] > 0 && used_node[get<2>(edge)] > 0) break;
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
      A[col]->insert_or_assign(k,d);
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

int main(int argc, char *argv[])
{
	// intT i;
	FILE *pFile=fopen(argv[1],"rb");
	const intT n=atol(argv[2]); //compile with -Wall flag
	int DIM=atoi(argv[3]);
  int rounds = atoi(argv[4]);
  string method = argv[5];
	float f;

	timer t;t.start();

	float *P = new float[n*DIM];
  double *D = new double[n*(n-1)/2];


	/************object generation*************/
	srand(time(0));
	for(intT i=0;i<n*DIM;i++)
	{
		if(1  == fscanf(pFile, "%f", &f)){
		P[i] = f;
		}else{
		printf("Failed to read point.\n");
		}
	}
	cout << "init " << t.next() << endl;
	cout << endl;

	for(int r = 0; r< rounds; ++ r){
		cout << "Round " << r << endl;
		timer t2;t2.start();
    parallel_for (intT i=0; i<n; i++) {
      parallel_for (intT j=i+1; j<n; j++) {
          double dist = 0;
          for(int d = 0; d < DIM; ++d){
            dist += (P[i*DIM + d] - P[j*DIM + d]) * (P[i*DIM + d] - P[j*DIM + d]);
          }
          if(method ==  "avgsq"){
            D[condensed_index(n,i,j)] = dist;
          }else{
            D[condensed_index(n,i,j)] = sqrt(dist);
          }
      // cout << i <<  " " << j << " " << condensed_index(n,i,j) << " " << sqrt(dist) << endl;
      }
    }

	cout << "distM " << t.next() << endl;

  sklearn(D, n, Complete());

  dendroentry *Z;
    if(method ==  "complete"){
      Z = nn_chain(D, n, Complete());
    }else if(method ==  "ward"){
      Z = nn_chain(D, n, Ward());
    }else if(method ==  "avg"){
      Z = nn_chain(D, n, Average());
    }else if(method ==  "avgsq"){
      Z = nn_chain(D, n, Average());
    }else{
        cout << "invalid method" << endl;
        exit(1);
    }
  
  free(Z);

	cout << "NNchain " << t.next() << endl;
	cout << "Total " << t2.next() << endl;
	cout << endl;

  }

	// for(intT i=0; i< 30; ++i){
	// 	cout << result[i] << endl;
  
	fclose(pFile);
	// t.reportTotal("total");
	return 0;
}
