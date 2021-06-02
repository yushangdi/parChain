#ifndef LINKAGE_DENDROFORMAT_H
#define LINKAGE_DENDROFORMAT_H

#include <fstream>
#include <iostream>

#include "sequence.h"
#include "sampleSort.h"
#include "gettime.h"


namespace dendrogram {

template<class nodeT>
struct nodeComparator{
  double eps;
  nodeComparator(double _eps):eps(_eps) {
  }

  bool operator() (nodeT &i, nodeT &j) {
      if(abs(i.getHeight() - j.getHeight()) <= eps) return i.getRound() < j.getRound();
     return i.getHeight() < j.getHeight();
  }
};

template<class nodeT>
struct DendroNode{
    intT cId;
    intT round = 0;
    intT idx; 
    intT n = 1;
    intT left = -1;
    intT right = -1;
    double height = 0;

    DendroNode(intT  t_cid, intT t_round, intT t_idx, intT _n, intT t_left, intT t_right, double _height):
        cId(t_cid),
        round(t_round),
        idx(t_idx),
        n(_n),
        left(t_left),
        right(t_right),
        height(_height){
    }
    
    DendroNode(nodeT &node):
        cId(node.cId),
        round(node.getRound()),
        idx(node.getIdx()),
        n(node.size()),
        left(node.left->getIdx()),
        right(node.right->getIdx()),
        height(node.getHeight()){
    }

    inline double getHeight(){return height;}
    inline intT getRound(){return round;}
    inline intT getIdx(){return idx;}
    inline intT size(){return n;}
};


struct dendroLine{
    intT id1;
    intT id2;
    double height;
    intT size;
    dendroLine(intT _id1, intT _id2, double _height, intT _size):id1(_id1), id2(_id2), height(_height), size(_size){}

    void print(ofstream file_obj){
        file_obj << id1 << " " << id2 << " " << std::setprecision(20) << height << " " << size << endl; 
    }

    void print(){
        cout << id1 << " " << id2 << " " << height << " " << size << endl; 
    }
};

template<class nodeT>
dendroLine* formatDendrogram(nodeT *nodes, intT n, double eps){
    DendroNode<nodeT> *sorted_nodes = newA(DendroNode<nodeT>, n-1);
    parallel_for(intT i=0; i<n-1; ++i){
        sorted_nodes[i] = DendroNode<nodeT>(nodes[i+n]);
    }

    sampleSort(sorted_nodes, n-1, nodeComparator<DendroNode<nodeT>>(eps));
    intT *map = newA(intT , n);
    parallel_for(intT i=0; i<n-1; ++i){
        map[sorted_nodes[i].getIdx() - n] = i+n;
    }
    dendroLine* dendrogram = newA(dendroLine, n-1);
    parallel_for(intT i=0; i<n-1; ++i){
        intT left = sorted_nodes[i].left;
        intT right = sorted_nodes[i].right;

        left = left < n ? left : map[left-n];
        right = right < n ? right : map[right-n];

        if(left > right) swap(left, right);


        dendrogram[i] = dendroLine(left, right,sorted_nodes[i].getHeight(),sorted_nodes[i].size());
    }

    free(map);
    free(sorted_nodes);
    return dendrogram;
}

}

#endif