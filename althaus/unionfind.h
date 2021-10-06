#ifndef UF_H
#define UF_H

#include <limits>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace parlay;

namespace UTIL {
  inline void PrintCaption(string t_in) {
    cout << "========= " << t_in << " =========" << endl;
  }

  inline void PrintSubcaption(string t_in) {
    cout << "========= " << t_in << endl;
  }

  template <class T>
  inline void PrintFunctionTitle(string t_func, T t_suffix) {
    cout << "[" << t_func << "] " << t_suffix << endl;
  }

  template <class T>
  inline void PrintFunctionItem(string t_func, string t_item, T t_suffix) {
    cout << "[" << t_func << "] " << t_item << " = " << t_suffix << endl;
  }

  template <class T>
  inline void PrintSubtimer(string t_item, T t_suffix, int prec = 5) {
    cout << "::" << t_item << ": " << std::setprecision(prec) << t_suffix << endl;
  }

  inline void PrintBreak() {
    cout << endl;
  }

  template <typename A>
  inline void PrintVec(A *t_vec, int t_len) {
      for (int i = 0; i < t_len; ++ i) {
	      cout << t_vec[i] << ' ';
      }
      cout << endl;
    }

  template <typename A>
  inline void PrintVec2(A *t_vec, int t_len) {
      for (int i = 0; i < t_len; ++ i) {
	     t_vec[i].print();
      }
      cout << endl;
    }

  template<class nodeT>
  inline void printNode(nodeT *Q, string s){
      cout << "===" << s << endl;
      (Q->pMin).print();
      (Q->pMax).print();
      cout << "===" << endl;
  }

  template<class nodeT>
  inline void printNode(nodeT *Q, int s){
      cout << "===" << s << endl;
      (Q->pMin).print();
      (Q->pMax).print();
  }

  template<class E>
  inline void printPair(E foo){
    std::cout << foo.first << ", " << foo.second << '\n';
  }

  template<class E>
  inline void printuple(E foo){
    std::cout << get<0>(foo) << ", " << get<1>(foo) << '\n';
  }

}

namespace UnionFind {
	
  template <typename IntType>
    struct ParUF {
			

    IntType *parents;
    pair<IntType, IntType> *hooks; // the edge that merged comp idx with a higher idx comp
    IntType m_n;
	double *values;
	bool store_value;
		
		ParUF() {}
      // initialize with all roots marked with -1
    ParUF(IntType n, bool t_store_value = false) {
		m_n = n;
		parents = (IntType *)malloc(sizeof(IntType)*n);
		parallel_for (0,n, [&](IntType i){parents[i] = std::numeric_limits<IntType>::max();}) ;

		hooks = (pair<IntType, IntType> *) malloc(sizeof(pair<IntType, IntType>) * n);
		parallel_for (0,n, [&](IntType i){hooks[i] = make_pair(std::numeric_limits<IntType>::max(), std::numeric_limits<IntType>::max());}) ;

		
		store_value = t_store_value;
		if (store_value){
			values = (double *)malloc(sizeof(double)*n);
		}else{
			values = nullptr;
		}
    }

      void del() {free(parents); free(hooks); if(values) free(values);}

      // Assumes root is negative 
      // Not making parent array volatile improves
      // performance and doesn't affect correctness
      inline IntType find(IntType i) {
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
			} 
			}
		}else{
			for(IntType i = 0; i < m_n; ++ i) {
			pair<IntType, IntType> ufEdge = get_edge(i);
			if (ufEdge.first < m_n) {
				file_obj << ufEdge.first << " ";
				file_obj <<  ufEdge.second << endl;
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

    };
  
} // end namespace UnionFind

#endif

