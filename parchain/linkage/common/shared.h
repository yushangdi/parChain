#ifndef SHARED_H
#define SHARED_H

#include <tuple>
#include "geometry.h"
// #include "graph.h"
#include "ndHash.h"
#include "serialHash.h"

#define LARGER_THAN_UB numeric_limits<double>::max()
#define CHECK_TOKEN -1
#define UNFOUND_TOKEN -2

// *************************************************************
//   Data Structure
// *************************************************************
template<int dim>
struct iPoint {
  typedef double floatT;
  typedef point<dim> pointT;
  intT i;
  pointT p;
  iPoint(pointT pp, intT ii): p(pp), i(ii) {}
  iPoint(pointT pp): p(pp), i(-1) {}
  iPoint(): i(-1) {}
  bool isEmpty() {return i<0;}
  floatT operator[](intT i) {return p[i];}
  floatT pointDist(iPoint q) {return p.pointDist(q.p);}
  floatT pointDist(iPoint *q) {return p.pointDist(q->p);}
  floatT dist(iPoint q) {return p.pointDist(q.p);}
  floatT pointDist(pointT q) {return p.pointDist(q);}
  floatT pointDistSq(iPoint q) {return p.pointDistSq(q.p);}
  floatT pointDistSq(iPoint *q) {return p.pointDistSq(q->p);}
  floatT pointDistSq(pointT q) {return p.pointDistSq(q);}
  intT idx() {return i;}
  void idx(intT ii) {i=ii;}
  pointT pt() {return p;}
  floatT* coordinate() {return p.x;}
  floatT coordinate(intT i) {return p.x[i];}
  floatT x(int i) {return p.x[i];}
  void x(int i, floatT val) {p.x[i]=val;}
  intT getDim(){return dim;}
  void minCoords(iPoint q){p.minCoords(q.p);}
  void maxCoords(iPoint q){p.maxCoords(q.p);}
  void print(){p.print();}
};

namespace LDS{// linkage data structure

enum Method { WARD, COMP, AVG, AVGSQ };

//linked list nodes
template <class eType>
struct node {
  eType elt;
  node* next;
node(eType _elt, node* _next): elt(_elt), next(_next) {}
};


struct EDGE{//nodes unordered
    volatile intT first;
    volatile intT second;
    volatile double w;//weight

    EDGE(intT t_u, intT t_v, double t_w):first(t_u), second(t_v), w(t_w){
        // if(first > second)swap(first,second);
    }
    EDGE():first(-1), second(-1), w(-1){}

    inline void print(){
        std::cout << "(" << first << ", " <<  second << "):" << w << std::endl;
    }

    inline double getW() const {return w;}
    inline pair<intT, intT> getE() const {return make_pair(first,second);}
    inline void update(intT t_u, intT t_v, double t_w){
      first = t_u;
      second = t_v;
      w = t_w;
    }
};

  // struct cInfoET {
  //   volatile intT idx;
  //   volatile double dist;
  //   cInfoET():idx(-1), dist(0){}
  //   cInfoET(intT ii, double dd):idx(ii), dist(dd){}

  //   bool operator== (const cInfoET& y) const
  //  {
  //      return idx == y.idx &&  dist == y.dist;

  //  }

  //  bool operator!= (const cInfoET& y) const
  //  {
  //      return idx != y.idx ||  dist != y.dist;

  //  }
  // };


  struct hashCluster { 
    typedef pair<intT,intT> eType;
    typedef intT kType;
    eType empty() {return pair<intT,intT>(-1,-1);}
    kType getKey(eType v) { return v.first; }
    uintT hash(intT s) { return utils::hash(s);}
    int cmp(intT v, intT b) {return (v > b) ? 1 : ((v == b) ? 0 : -1);}
    bool replaceQ(eType s, eType s2) {cout << "wrong replace Q" << endl;
        exit(1);return 0;}//return s.second > s2.second;}
    bool replaceQ(eType s, eType* TA, intT h) {
        utils::writeAdd(&(TA[h].second), s.second);
    return false;}
  };

  struct hashClusterAveET {
    volatile intT first;
    volatile intT idx;
    volatile double dist;
    hashClusterAveET():first(-1), idx(-1), dist(UNFOUND_TOKEN){}
    // hashClusterAveET(intT ii, cInfoET jj):first(ii), second(jj){}
    hashClusterAveET(intT ii, intT jj, double dd):first(ii), idx(jj), dist(dd){}
    void print(){
      cout << first << " " << idx << " " << dist << endl;
    }

    bool operator== (const hashClusterAveET& y) const
   {    
       return first == y.first && idx == y.idx &&  abs(dist - y.dist) < 1e-20; //todo: CHANGE? should only be used for detect empty

   }

       bool operator!= (const hashClusterAveET& y) const
   {   
       return first != y.first || idx != y.idx ||  abs(dist - y.dist) > 1e-20;

   }

  };

  struct hashClusterAve { 
      typedef hashClusterAveET eType; // (cid, (idx, dist))
      typedef intT kType;
      eType empty() {return eType(-1,-1,UNFOUND_TOKEN);}
      kType getKey(eType v) { return v.first; }
      intT hash(intT s) { return utils::hashInt(s);}
      int cmp(intT v, intT b) {return (v > b) ? 1 : ((v == b) ? 0 : -1);}
      bool replaceQ(eType s, eType s2) { // s is the new value
        if(s.dist == CHECK_TOKEN){ // insert check
          if(s.idx == s2.idx){ return false;} // already exist, no need to compute
          else{return true;}// wrong index, need to compute
        }
        // a real insert, replace
        return true;
        // if non TOKEN s is updating in rangeQuery, current entry must have been inserted TOKEN
        // or already exist  with valid idx, so this update cannot be overwritten by a TOKEN update
      }//return s.second > s2.second;}
      bool replaceQ(eType s, eType* TA, intT h, eType s2) {
        if(s.dist == CHECK_TOKEN){ // insert check
          if(s.idx == s2.idx){ return false;} // already exist, no need to compute
          else{return true;}// wrong index, need to compute
        }
        // a real insert, replace
        return true;
      }
  };

using distCacheT = Table<hashClusterAve, intT>;
using clusterCacheT = SerialTable<hashCluster,intT >;

struct edgeComparator{
    bool operator () (EDGE i, EDGE j) {
        return i.getW() < j.getW();
        }
};

struct edgeComparator2Debug{
    double eps = 1e-20;
    edgeComparator2Debug(){}
    edgeComparator2Debug(double _eps):eps(_eps){}
    bool operator () (EDGE i, EDGE j) {
      // if(abs(i.getW() - j.getW()) < i.getW()*std::numeric_limits<double>::epsilon()) return i.second < j.second;
      if(i.first == 41705){
        cout << "=======" << endl;
        cout << i.second <<   " " << j.second << endl;
        cout << std::setprecision(23) << i.getW() << endl;
        cout << std::setprecision(23) << j.getW() << endl;
        cout <<  std::setprecision(23) << i.getW() - j.getW() << endl;
        cout << (abs(i.getW() - j.getW()) <= eps) << endl;
        cout << "=======" << endl;
      }
      if(abs(i.getW() - j.getW()) <= eps) return i.second < j.second;
      return i.getW() < j.getW();
      }
};

struct edgeComparator2{
    double eps = 1e-20;
    edgeComparator2(){}
    edgeComparator2(double _eps):eps(_eps){}
    bool operator () (EDGE i, EDGE j) {
      // if(abs(i.getW() - j.getW()) < i.getW()*std::numeric_limits<double>::epsilon()) return i.second < j.second;
      if(abs(i.getW() - j.getW()) <= eps) return i.second < j.second;
      return i.getW() < j.getW();
      }
};

struct edgeComparatorMax{
    bool operator () (EDGE i, EDGE j) {
        return (i.getW() > j.getW());
        }
};

template<class PT>
struct PairComparator2{
  bool operator () (PT i, PT j){
    return i.second < j.second;
  }
};

template<class PT>
struct PairComparator21{
  double eps = 1e-20;
  PairComparator21(){}
  PairComparator21(double _eps):eps(_eps){}
  bool operator () (PT i, PT j){
    // return (i.second < j.second) || ((i.second ==  j.second) && (i.first < j.first));
    // if(abs(i.second - j.second) < i.second*std::numeric_limits<double>::epsilon()) return i.first < j.first;
    if(abs(i.second - j.second) <= eps) return i.first < j.first;
      return i.second < j.second;
  }
};

// sample sort struct
struct ClusterIdComparator {
  intT *m_clusterIds;
  ClusterIdComparator(intT *t_clusterIds) {
    m_clusterIds = t_clusterIds;
  }
  bool operator() (const intT i, const intT j) {
    return m_clusterIds[i] < m_clusterIds[j];
  }
};

template<class pointT>
inline double pointDist(pointT* p, pointT* q){
  return p->pointDist(q);
}

template<class pointT>
inline double pointDist(pointT p, pointT q){
  return p.pointDist(q);
}

template<class pointT>
inline double getPointId(pointT* p){
  return p->idx();
}

template<class pointT>
inline double getPointId(pointT p){
  return p.idx();
}

}

// *************************************************************
//   Misc
// *************************************************************

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
      for (intT i = 0; i < t_len; ++ i) {
	      cout << t_vec[i] << ' ';
      }
      cout << endl;
    }

  template <typename A>
  inline void PrintVec2(A *t_vec, int t_len) {
      for (intT i = 0; i < t_len; ++ i) {
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
  inline void printNode(nodeT *Q, intT s){
      cout << "===" << s << endl;
      (Q->pMin).print();
      (Q->pMax).print();
  }

  template<class E>
  inline void printPair(E foo){
    std::cout << foo.first << ", " << foo.second << '\n';
  }

  template<class E>
  inline void printTuple(E foo){
    std::cout << get<0>(foo) << ", " << get<1>(foo) << '\n';
  }

  template<class tableT>
  inline void PrintTable(tableT *tb){
      for(intT i = 0; i < tb->m; ++i){
        if(tb->TA[i] != tb->empty){
          (tb->TA[i]).print();
        }
      }
  }

  template<int dim, class nodeT>
  struct  PrintClusterId{

      PrintClusterId(){}

      inline bool isTopDown(intT id){ return true;}

      inline void TopDownNode(nodeT *Q, intT id){
          UTIL::printNode(Q, Q->nInfo.getCId());
      }

      inline void BottomUpNode(nodeT *Q, intT id){}

      inline void BaseCase(nodeT *Q, intT id){
          UTIL::printNode(Q, Q->nInfo.getCId());
      }

      inline intT SwitchMode(nodeT *Q, intT id){return -1;}

      inline bool Par(nodeT *Q){return false;}

      inline bool Stop(nodeT *Q, intT id){return false;}

  };



}

#endif