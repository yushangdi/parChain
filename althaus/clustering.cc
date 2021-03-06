#include<iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <set>
#include <queue>

#include "gettime.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/range.h"

#include "parlay/slice.h"
#include "parlay/primitives.h"
#include "unionfind.h"


using namespace std;	
using namespace parlay;

using T=float;

// #define MAX_T
#define MAX_T numeric_limits<T>::max()

class pointset {
  private:
    // we assume points in R^dimension
    int dimension;

    // number of points in Pointset
    int size_;

    // Points[i][j] gives the jth coordinate of pint i 
    T** Points;

  public:

    pointset(int dimension_) {
      dimension=dimension_;
      size_=0;
    }

    void from_file(FILE *pFile, int s) {
      size_ = s;
      T h;
      Points=new T*[size_];
      for(int i=0; i<size_; i++) {
        Points[i]=new T[dimension];
        for(int j=0; j<dimension; j++) {
          if(1  == fscanf(pFile, "%f", &h)){
          Points[i][j]=h;
          }else{
          printf("Failed to read point.\n");
          }
        }
      }
    }

    // we create a pointset of the given size of point sampled uniformely at random
    // We assume that this function is called exactly once
    void uniform_random(int s) {
      size_=s;
      Points=new T*[size_];
      for(int i=0; i<size_; i++) {
       Points[i]=new T[dimension];
       for(int j=0; j<dimension; j++) {
         Points[i][j]=(T) random()/RAND_MAX;
       }
      }
    }

    // We compute the distance of the points i and j
    T point_dist(int i, int j) const {
      T d=0;
      for(int k=0; k<dimension; k++) d+=(Points[i][k]-Points[j][k])*(Points[i][k]-Points[j][k]);
      return sqrt(d);
    }

    // Return number of Points in the set
    int size() const {
      return size_;
    }

    // we delete the arrays in which the points were stored
    // We do not check (but assume), whether there were created  
    ~pointset() {
      for(int i=0; i<size(); i++) delete[] Points[i];
      delete[] Points;
    }
};

    // We compute the distance of the points i and j
    template<class PS>
    T point_dist(const PS& P, int i, int j) {
      return P.point_dist(i,j);
    }
 

 class myPQ{

   public:
   typedef tuple<T, int> entry ;
   vector<entry> v;

  void pop(){
    pop_heap(v.begin(), v.end());
    v.pop_back();
  }

    entry top(){
      entry tmp = v.front();
      get<0>(tmp) = -1 * get<0>(tmp);
      return tmp;
    }

    void push(entry i){
      entry tmp = i;
      get<0>(tmp) = -1 * get<0>(tmp);
      v.push_back(tmp);
      
    // using push_heap() to reorder elements
    push_heap(v.begin(), v.end());
    }

    void push_to_vec(entry i){
      entry tmp = i;
      get<0>(tmp) = -1 * get<0>(tmp);
      v.push_back(tmp);
      
    }

    auto begin(){
      return v.begin();
    }
    auto end(){
      return v.end();
    }

    auto size(){return v.size();}


 };


// use a set to store the removed element, check upon pop
class junk_array {
  typedef tuple<T, int> entry ;
  typedef myPQ pq;
  
  private:

    int size;
    int num_processor;
    int num_distances;

    UnionFind::ParUF<int> uf;
    set<int> active_clusters;

    sequence<entry> help;
    pq **queues;
    int* Mapping;//any point id of a point in the cluster, or -1 if the cluster is merged. cluster id to uf id
    int* Mapping2; // uf id to clsuter id

  public:

    junk_array(pointset &P, int num_distances_) {
      size=P.size();
      num_processor = num_workers();
      num_distances=num_distances_;
      uf = UnionFind::ParUF<int>(size);

      help = sequence<entry>(size * num_processor);

      Mapping=new int[2*size];
      Mapping2=new int[size];

      for(int i=0; i<size; i++) Mapping[i]=i;
      for(int i=size; i<2*size; i++) Mapping[i]=-1;
      for(int i=0; i<size; i++) Mapping2[i]=i;

      
      queues = (pq **)malloc(sizeof(pq *) * 2* size);

      parallel_for (0,size, [&](int i){
        queues[i] = new pq();
        fillQueue(P, i);
      }) ;

      for(int i=0; i<size; i++){
        active_clusters.insert(i);
      }

    }

    ~junk_array() {
      delete[] Mapping;
      delete[] Mapping2;

      delete[] queues;
    }

    int current_cluster(int i) {
      if(i >= size) return uf.find(Mapping[i]);
      return uf.find(i);
    }

    void fillQueue(pointset &P, int i)
    {
      
      int help_start = worker_id() * size;
      for(int j=0; j<P.size(); j++) { help[j + help_start]=tuple<T, int>(-1,-1); }
      
      for(int ii=0; ii<P.size(); ii++) if(current_cluster(ii)==current_cluster(i)) { // loop over points in the same cluster as i 
        for(int j=0; j<P.size(); j++) { // loop over other clusters
          int k=current_cluster(j);
          if(k == current_cluster(i)) continue;
          T d=point_dist(P,ii,j);
          if(d>get<0>(help[k+ help_start])) { //only store one copy does not work, need both copies when merging
            if(k>P.size()) std::cout<<"ERROR 2\n";
            help[k+ help_start]=tuple<T, int>(d, Mapping2[k]);
          }
        }
      }

      auto help_valid = parlay::filter(make_slice(help).cut(help_start, size + help_start), [&](entry i){return get<1>(i)!= -1;});

      if(help_valid.size()== 0) return;
      // Sort the given array
      parlay::sort_inplace(help_valid);

      for (int ii=0;ii<min(num_distances, (int)help_valid.size());++ii){
        if(get<1>(help_valid[ii])!= -1) queues[i]->push_to_vec(help_valid[ii]);
      }
      if(queues[i]->size()> 0) make_heap(queues[i]->begin(),queues[i]->end());
    }

    bool no_distances(int i) {
      return queues[i]->size() == 0;
    }

    void pop(pq *v){
      v->pop();
    }

    int peek_val(pq *v){
      return get<1>(v->v.front());

    }
    T peek_key(pq *v){
      return -1 * get<0>(v->v.front());
    } 

    // Compute the mimial distance of a pair of clusters and return the clusters in i and j
    // Todo so, we iterate over all junks. If the junk is used (i.e. has at least one distance stored),
    // we update the smallest distance found so far
    T get_minimal_distance(pointset &P, int& i, int& j) {
      // std::cout<<"look for miminal distance\n";
      T d=MAX_T;
      for(auto itr = active_clusters.begin(); itr != active_clusters.end(); itr++){
        int k = *itr;
        pq *v = queues[k];
        while(!no_distances(k) &&  (peek_val(v)== -1 || Mapping[peek_val(v)] == -1)){ //
          pop(v); // skip invalid entries
        }
        if(no_distances(k)) {
          fillQueue(P, k);
        }
        if(!no_distances(k) && v->top() < make_tuple(d,j)){
          d = peek_key(v);
          i = k;
          j = peek_val(v);
        }
      }
      pop(queues[i]);
      return d;
    }

    // remove the smaller of j and k from queues[i]
    void remove(int j, int k, int i, int newid){
      // if only one of them in the queue, remove
      //if both in the queue, remove one, and update the val of the other to be newid
      pq *v = queues[i];
      vector<entry>::iterator first_itr = v->end();
      for(auto it = v->begin(); it != v->end(); ++it) {
        
        if(get<1>(*it) == j || get<1>(*it) == k){
          if(first_itr == v->end()){
            first_itr = it;
          }else{
            if(get<0>(*it) > get<0>(*first_itr)){
              get<1>(*first_itr)=newid; //  first is larger
              get<1>(*it)=-1;
            }else{
              get<1>(*it)=newid; //this is larger
              get<1>(*first_itr)=-1;
            }
            first_itr = v->end();
            break;
          }
        }
      }

       if(first_itr != v->end()){
         get<1>(*first_itr)=-1;
       }

    }

    void merge(int k, int j, int i) {
      queues[k] = new pq();
      pq *vi = queues[i];
      pq *vj = queues[j];

      auto min1 = min_element(vi->v);
      auto min2 = min_element(vj->v);
      T upper = -1 * max(get<0>(*min1), get<0>(*min2));

      for(auto it = vi->begin(); it != vi->end(); ++it) {
        int key = get<1>(*it);
        if(key == -1 || Mapping[key]== -1) continue;
       for(auto it2 = vj->begin(); it2 != vj->end(); ++it2) {
        if(key == get<1>(*it2)){
          T newd = -1 * min(get<0>(*it), get<0>(*it2));
          if(newd < upper) queues[k]->push_to_vec(entry(newd, key)); 
          break;
        }
        }
      }
      if(queues[k]->size()> 0) make_heap(queues[k]->begin(),queues[k]->end());
      // if we want to update other queues, only need to check existing ones, and update, because if one is not in the queue, it can only be larger
      delete queues[i];
      delete queues[j];
      uf.link(Mapping[i],Mapping[j]);
      Mapping[k] = uf.find(Mapping[i]);
      Mapping2[uf.find(Mapping[i])] = k;
      Mapping[i] = -1;
      Mapping[j] = -1;
      active_clusters.erase(i);
      active_clusters.erase(j);
      active_clusters.insert(k);

      int* tmp_active = new int[active_clusters.size()];
      int tmp_i=0;
      for(auto itr = active_clusters.begin(); itr != active_clusters.end(); itr++){
        int r = *itr;
        tmp_active[tmp_i] = r;
        tmp_i++;
      }
      parallel_for(0, (int)active_clusters.size(), [&](int itr){
        int r = tmp_active[itr];
        remove(i,j,r, k);
      }) ;
      delete[] tmp_active;

    };

};


int main(int argc, char *argv[]) {
  FILE *pFile=fopen(argv[1],"rb");
	const int MAX_N=atol(argv[2]); 
	int DIM=atoi(argv[3]);
  int num_dist = atoi(argv[4]);

  pointset P(DIM);
  P.from_file(pFile, MAX_N);
  fclose(pFile);
  timer t; t.start();
  std::cout<<"Clustering "<<P.size()<<" points\n";
  cout << "using " <<num_dist << " cache" << endl;
  T cost = 0;

  cout << "using " << parlay::num_workers() << " threads" << endl;

  junk_array JA = junk_array(P, num_dist);

  for(int i=P.size(); i<2*P.size()-1; i++) {
    int j,k;
    T d=JA.get_minimal_distance(P,j,k);
    cost += d;
    if(j==k) { std::cout<<"ERROR i==j\n"; exit(1); }
    JA.merge(i,j,k);
  }

  cout << "cost: " << std::setprecision(8) << cost << endl;
  cout << "time: " << t.next() << endl;

  return 0;
};



