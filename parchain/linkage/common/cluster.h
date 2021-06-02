#ifndef CLUSTER_H
#define CLUSTER_H

#include "unionfind.h"
#include "sequence.h"
#include "sampleSort.h"

// sample sort struct
template<class objT>
struct ClusterIdComparator {
  UnionFind::ParUF<intT> *uf;
  ClusterIdComparator(UnionFind::ParUF<intT> *t_uf) {
    uf = t_uf;
  }

  bool operator() (objT &i, objT &j) {
    return uf->find(i->idx()) < uf->find(j->idx());
  }
};

template<class objT>
inline pair<intT, intT*> GetCluster_helper(UnionFind::ParUF<intT> *t_uf, objT **t_out, intT t_n) {

  sampleSort(t_out, t_n, ClusterIdComparator<objT *>(t_uf));

  intT *flags = newA(intT, t_n);
  flags[0] = 1;
  parallel_for (intT i = 1; i < t_n; ++ i) {
    if (t_uf->find(t_out[i]->idx()) != t_uf->find(t_out[i - 1]->idx())) {
      flags[i] = 1;
    }else{
      flags[i] = 0;
    }
  }

  intT numClusters;
//   if (sizeof(intT) > 4) {
//     numClusters = sequence::prefixSumL<intT>(flags, 0, t_n);
//   } else {
    numClusters = sequence::prefixSum<intT>(flags, 0, t_n);
//   }
  intT *clusterOffsets = newA(intT, numClusters + 1);
  parallel_for (intT i = 0; i < t_n - 1; ++ i) {
    if (flags[i] != flags[i + 1]) {
      clusterOffsets[flags[i]] = i;
    }
  }
  intT i = t_n - 1;
  if (flags[i] != numClusters) {
    clusterOffsets[flags[i]] = i;
  }
  clusterOffsets[numClusters] = t_n;

  free(flags);
  return pair<intT, intT*>(numClusters, clusterOffsets);
}

template<class objT>
inline pair<intT, intT*> GetCluster(UnionFind::ParUF<intT> *t_uf, objT **t_out, intT *t_activeClusters, intT t_n) {
  pair<intT, intT *> myClusters = GetCluster_helper<objT>(t_uf, t_out, t_n);
  parallel_for(intT i = 0; i < myClusters.first; ++ i){
      t_activeClusters[i] = t_uf->find(t_out[myClusters.second[i]]->idx());
  }
  return myClusters;
}



// template <typename intT>
// class ActivePointFilter
// {
// public:
//   intT *PMapping;
//   double epsilon;
//   double * bcpDists;
//   ActivePointFilter(intT *t_mapping, double t_epsilon, double *t_bcpDists) :
//     PMapping(t_mapping), epsilon(t_epsilon), bcpDists(t_bcpDists) {};
//   bool operator()(intT i)
//   {
//     return bcpDists[PMapping[i]] <= epsilon;
//   }
// };


// //(myUF, PIndices, PMapping, n);
// template<class intT>
// inline pair<intT, intT *> GetActiveCluster(UnionFind::ParUF<intT> *t_uf, intT *t_out, intT *t_mapping, intT t_n, intT &n_active, double epsilon, double *bcpDists) {

//     intT new_n_active = distance(t_out, partition(t_out, t_out + n_active, ActivePointFilter<intT>(t_mapping, epsilon, bcpDists)));
//     n_active = new_n_active;

//     intT *clusterIds = newA(intT, t_n);
//     parallel_for(intT i = 0; i < n_active; ++ i) {
//       clusterIds[t_out[i]] = t_uf->find(t_out[i]);
//     }

//     pair<intT, intT *> myClusters = GetCluster_helper<intT>(t_uf, t_out, t_mapping, n_active, clusterIds);

//     free(clusterIds);
//     return myClusters;
// }

// template<class intT, class F>
// inline pair<intT, intT *> GetActiveCluster(UnionFind::ParUF<intT> *t_uf, intT *t_out, intT *t_mapping, intT t_n, intT &n_active, double epsilon, F f) {

//     intT new_n_active = distance(t_out, partition(t_out, t_out + n_active, f));
//     n_active = new_n_active;

//     intT *clusterIds = newA(intT, t_n);
//     parallel_for(intT i = 0; i < n_active; ++ i) {
//       clusterIds[t_out[i]] = t_uf->find(t_out[i]);
//     }

//     pair<intT, intT *> myClusters = GetCluster_helper<intT>(t_uf, t_out, t_mapping, n_active, clusterIds);

//     free(clusterIds);
//     return myClusters;
// }

#endif
