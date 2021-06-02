#ifndef UF_H
#define UF_H

#include <limits>
#include <fstream>
#include <iostream>
#include <string>
#include "shared.h"

using namespace std;

namespace UnionFind {
  template <typename IntType>
    
    struct ParUF {

    IntType *parents;
    pair<IntType, IntType> *hooks; // the edge that merged comp idx with a higher idx comp
    IntType m_n;
	double *values;
	bool store_value;

      // initialize with all roots marked with -1
    ParUF(IntType n, bool t_store_value = false) {
		m_n = n;
		parents = newA(IntType, n);
		parallel_for (IntType i=0; i < n; i++) parents[i] = std::numeric_limits<IntType>::max();
		hooks = (pair<IntType, IntType> *) malloc(sizeof(pair<IntType, IntType>) * n);
		parallel_for (IntType i=0; i < n; i++) {
		hooks[i] = make_pair(std::numeric_limits<IntType>::max(), std::numeric_limits<IntType>::max());
		}
		
		store_value = t_store_value;
		if (store_value){
			values = newA(double, n);
		}else{
			values = nullptr;
		}
    }

      void del() {free(parents); free(hooks); if(values) free(values);}

      // Assumes root is negative 
      // Not making parent array volatile improves
      // performance and doesn't affect correctness
      inline IntType find(IntType i) {
				// if(i >= m_n){
				// 	cout << "uf find larger than init size" << endl;
				// 	return i;
				// }
	IntType j = i;
	if (parents[j] == std::numeric_limits<IntType>::max()) return j;
	do j = parents[j];
	while (parents[j] < std::numeric_limits<IntType>::max());
	//note: path compression can happen in parallel in the same tree, so
	//only link from smaller to larger to avoid cycles
	IntType tmp;
	while((tmp=parents[i])<j){ parents[i]=j; i=tmp;} 
	return j;
      }

    IntType link(IntType u, IntType v) {
		IntType c_from = u;
		IntType c_to = v;
		while(1){
		u = find(u);
		v = find(v);
		if(u == v) break;
		if(u > v) swap(u,v);
		//if successful, store the ID of the edge used in hooks[u]
		if(hooks[u].first == std::numeric_limits<IntType>::max() && __sync_bool_compare_and_swap(&hooks[u].first,std::numeric_limits<IntType>::max(),c_from)){
			parents[u]=v;
			hooks[u].second=c_to;
			break;
		}
		}
		return parents[u];
    }

	IntType link(IntType u, IntType v, double lv) {
		IntType c_from = u;
		IntType c_to = v;
		while(1){
		u = find(u);
		v = find(v);
		if(u == v) break;
		if(u > v) swap(u,v);
		//if successful, store the ID of the edge used in hooks[u]
		if(hooks[u].first == std::numeric_limits<IntType>::max() && __sync_bool_compare_and_swap(&hooks[u].first,std::numeric_limits<IntType>::max(),c_from)){
			parents[u]=v;
			hooks[u].second=c_to;
			values[u] = lv;
			break;
		}
		}
		return parents[u];
    }

    pair<IntType, IntType> get_edge(IntType idx) {
		return hooks[idx];
    }

	pair<pair<IntType, IntType>, double> get_edge_value(IntType idx) {
		return make_pair(hooks[idx], values[idx]);
    }

	inline void print_edges(){
		if(store_value){
			for(IntType i = 0; i < m_n; ++ i) {
			pair<IntType, IntType> ufEdge = get_edge(i);
			if (ufEdge.first < m_n) {
				cout << ufEdge.first << " ";
				cout <<  ufEdge.second << " ";
				if(values) cout << std::setprecision(6)  <<  values[i] << endl;
				// cout << "========= " << endl;
			} 
			}
		}else{
			for(IntType i = 0; i < m_n; ++ i) {
			pair<IntType, IntType> ufEdge = get_edge(i);
			if (ufEdge.first < m_n) {
				cout << ufEdge.first << " ";
				cout <<  ufEdge.second << endl;
				cout << "========= " << endl;
			} 
			}
		}

	}

	inline void print_edges(string filename){
		ofstream file_obj;
		file_obj.open(filename);
		if(store_value){
			for(IntType i = 0; i < m_n; ++ i) {
			pair<IntType, IntType> ufEdge = get_edge(i);
			if (ufEdge.first < m_n) {
				file_obj << ufEdge.first << " ";
				file_obj <<  ufEdge.second << " ";
				if(values) file_obj << std::setprecision(6)  <<  values[i] << endl;
				// cout << "========= " << endl;
			} 
			}
		}else{
			for(IntType i = 0; i < m_n; ++ i) {
			pair<IntType, IntType> ufEdge = get_edge(i);
			if (ufEdge.first < m_n) {
				file_obj << ufEdge.first << " ";
				file_obj <<  ufEdge.second << endl;
				// cout << "========= " << endl;
			} 
			}
		}

	}

	inline double cost(){
		double result = 0;
		IntType ct = 0;
		for(IntType i = 0; i < m_n; ++ i) {
			pair<IntType, IntType> ufEdge = get_edge(i);
			if (ufEdge.first < m_n) {result += values[i]; ct += 1;}  
		}
		UTIL::PrintFunctionItem("CLINK", "Edge Num", ct );
		return result;
	}

	intT *fcluster(double eps, bool output = false, string fileName = "./cluster.txt"){

     	 UnionFind::ParUF<intT> *uf2 = new UnionFind::ParUF<intT>(m_n, false);
	 	 intT *clusters = newA(intT, m_n);

	  	parallel_for(IntType i = 0; i < m_n; ++ i) {
		  if(get_edge(i).first < m_n && values[i] <= eps){
			uf2->link(get_edge(i).first, get_edge(i).second);
		  }
		}
		
		parallel_for(intT i = 0; i < m_n; ++i){
      		clusters[i] = uf2->find(i);
    	}

		if(output){
			ofstream file_obj;
			file_obj.open(fileName);
			for (IntType i = 0; i < m_n; ++ i){
			file_obj << clusters[i] <<  ", ";
			}
			file_obj << endl;
		}

		delete uf2;

		return clusters;

	}	

	void serialize(string fileName = "./uf.txt"){
		ofstream file_obj;
		file_obj.open(fileName);
		file_obj << m_n << endl;
		file_obj << store_value << endl;
		for(intT i=0; i<m_n; ++i){
			file_obj << hooks[i].first << endl;
		}
		for(intT i=0; i<m_n; ++i){
			file_obj << hooks[i].second << endl;
		}
		for(intT i=0; i<m_n; ++i){
			file_obj << parents[i] << endl;
		}
		if(!store_value) return;
		for(intT i=0; i<m_n; ++i){
			file_obj << values[i] << endl;
		}
		file_obj.close();
	}

	void deserialize(string fileName = "./uf.txt"){
		ifstream file_obj;
		file_obj.open(fileName.c_str(), ios::in);
		string ns;
		getline(file_obj, ns);
		assert(m_n <= stoi(ns)); // otherwise will overflow

		getline(file_obj, ns, '\n');
		bool tmp = (stoi(ns) != 0);
		//assert(store_value == tmp); // otherwise will overflow

		for(intT i=0; i<m_n; ++i){
			getline(file_obj, ns, '\n');
			hooks[i].first =  stoi(ns);
		}
		for(intT i=0; i<m_n; ++i){
			getline(file_obj, ns, '\n');
			hooks[i].second =  stoi(ns);
		}
		for(intT i=0; i<m_n; ++i){
			getline(file_obj, ns, '\n');
			parents[i] =  stoi(ns);
		}
		if(store_value){
		for(intT i=0; i<m_n; ++i){
			getline(file_obj, ns, '\n');
			values[i] =  stod(ns);
		}
		}

		file_obj.close(); //close the file object.

	}
	

    };
  
} // end namespace UnionFind

#endif

