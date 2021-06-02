#pragma once

#include "sequence.h"
#include "kdTree2.h"
#include "shared.h"

using namespace std;

// *************************************************************
//   NN Chain
// *************************************************************
// remove chainRev for ave?
template<int dim>
struct TreeChainInfo{
  intT *terminal_nodes;// must be cluster id
  intT chainNum;
  intT *chain;//chain[i] is cluster i's nn, -1 for unknown, -2 for invalid
  pair<long, double> *chainRev;// reverse chain, indexed by cid, TODO: change to LDS
  bool *findNN;
  bool *flag;//used in next() to update terminal_nodes
  pair<intT, double> *revTerminals; // indexed by ith temrinal nodes
  LDS::PairComparator21<pair<long, double>> EC2;
  double eps;
  // const LDS::PairComparator2<pair<long, double>> EC2 = LDS::PairComparator2<pair<long, double>>();

  TreeChainInfo(intT n, double eps){
    EC2 = LDS::PairComparator21<pair<long, double>>(eps);
    terminal_nodes = newA(intT, n);
    chain = newA(intT, n);
    chainRev = (pair<long, double> *)malloc(sizeof(pair<long, double>) * n);
    findNN = newA(bool, n);
    flag = newA(bool, n);
    chainNum = n;
    revTerminals = (pair<intT, double> *)malloc(sizeof(pair<intT, double>) * chainNum);

    parallel_for(intT i = 0; i < n; ++ i){
      terminal_nodes[i] = i;
      chain[i] = -1;
      findNN[i] = true;
      chainRev[i] = make_pair((long)-1,numeric_limits<double>::max());
    }
  }

  ~TreeChainInfo(){
    free(terminal_nodes);
    free(chain);
    free(findNN);
    free(chainRev);
    free(flag);
    free(revTerminals);
  }


  inline void updateChain(intT cid, intT nn, double w){
    chain[cid] = nn;
    utils::writeMin(&chainRev[nn], make_pair((long)cid,w), EC2); 
  }

  //TODO: change chain[prev[cid]] to -1 as well if code = -1
  // then in checking, only check for -1
  inline void invalidate(intT cid, intT code){
    chain[cid] = code;
  }

  inline void invalidateRev(intT cid){
    chainRev[cid] = make_pair((long)-1,numeric_limits<double>::max());
  }

  //cid >= 0
  inline intT getNN(intT cid){
    return chain[cid];
  }
  //get the rev of ith terminal nodes 
  inline pair<intT, double> getChainPrev(intT i){
    return revTerminals[i];
  }

  // update findNN, terminal_nodes and chainNum
  template<class F>
  inline void next(F *finder){
#ifdef DEBUG
  UTIL::PrintVec(chain, finder->n);
#endif
    parallel_for(intT i = 0; i < finder->n; ++i){
      findNN[i] = false;
    }
    intT C = finder->C;
    parallel_for(intT i = 0; i < C; ++i){
      intT cid = finder->activeClusters[i];
      flag[i] = getNN(cid) == -1  || getNN(getNN(cid)) < 0;// ok because -2 won't be in active clusters
      findNN[cid] = flag[i];
      //become a terminal node if it's a new merged cluster or its nn is merged
    }
    _seq<intT> newClusters = sequence::pack<intT, intT, sequence::getA<intT,intT> >(NULL, flag, 0, C, sequence::getA<intT,intT>(finder->activeClusters));
    free(terminal_nodes);
    terminal_nodes = newClusters.A;
    chainNum = newClusters.n;

    parallel_for(intT i=0; i<chainNum; ++i){
      revTerminals[i] =  make_pair((intT)chainRev[terminal_nodes[i]].first, chainRev[terminal_nodes[i]].second);
    }
  }

};

template<int dim, class TF>
inline void check_chain(TF *finder){
    for(intT i=0;i<finder->C;++i){
      intT cid = finder->activeClusters[i];
      double w1 = finder->edges[cid].getW();
      intT cid2 = finder->edges[cid].second;
      double w2 = finder->edges[cid2].getW();
      if(w2 - w1 > 0){ //finder->eps
        cout << cid << " " << finder->edges[cid].second << " " << std::setprecision(20) << w1 << endl;
        cout << cid2 << " " << finder->edges[cid2].second << " " << std::setprecision(20) << w2 << endl;
        cout << "========" << endl;
      }
    }
}

template<class TF, class CH>
inline void zero_chain_debug(TF *finder, int round, CH* info){
    cout << "0 chain num !!!!!" << endl;
    intT cid = finder->activeClusters[0];
      for(intT i=0;i<100;++i){
        // cout << info->chain[cid] << endl;
        cout << cid << " " << finder->edges[cid].second << " " << std::setprecision(23) << finder->edges[cid].getW() << endl;
        cid = info->chain[cid];
    }
    exit(1);
}

// checking if chain is pointing to nearest neighbor
template<class TF, class CH>
inline void chain_debug(TF *finder, int round, CH* info){
    
      for(intT i=0;i<finder->C;++i){
        intT cid = finder->activeClusters[i];
        
        if(info->chain[cid] >= 0 && info->chain[cid] != finder->edges[cid].second){
        cout << cid << " " << info->chain[cid] << " " << finder->edges[cid].second << endl;
         exit(1);
         }
    }
   
}

template<int dim>
inline intT count_chain_length(TreeChainInfo<dim> *info, intT cid){
    intT nn = -1;
    intT ct = 0;
    while(cid != nn){
        ct++;
        cid = info->chain[cid];
        if(cid < 0) break;
        nn = info->chain[cid];
    }
    return ct;
}