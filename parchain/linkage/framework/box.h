#ifndef PACCHAIN_BOX_H
#define PACCHAIN_BOX_H

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include "shared.h"

namespace FINDNN {
  
  template<intT dim, class nodeT>
  struct queryBallSimple{
      inline double getBall(nodeT* query, double r){
        return r;
      }
  };

  template<intT dim, class nodeT>
  struct queryBallSqrt{
      inline double getBall(nodeT* query, double r){
        return sqrt(r);
      }
  };

  template<intT dim, class nodeT>
  struct queryBallWard{

      inline double getBall(nodeT* query, double r){
        double n = (double) query->size();
        return r*sqrt((n+1)/2.0/n);
      }
  };

}
#endif