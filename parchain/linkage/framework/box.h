#ifndef PACCHAIN_BOX_H
#define PACCHAIN_BOX_H

#include <limits>
#include <fstream>
#include <string>
#include <vector>
#include "shared.h"

namespace FINDNN {
  // find points in `root` within box of side `r` centered at `query` point
  // template<intT dim, class nodeT>
  // struct queryBoxComp1{

  //     inline tuple<point<dim>, point<dim>> getBox(nodeT* query, double r){
  //       point<dim> pMin1 = point<dim>();
  //       point<dim> pMax1 = point<dim>();
  //       // double* center = query->items[0]->coordinate();
  //       for (int i=0; i<dim; ++i) {
  //       double center_i = (query->pMin[i] + query->pMax[i])/2; // distance to center
  //       pMin1.updateX(i, center_i - r);
  //       pMax1.updateX(i, center_i + r);}
  //       return make_tuple(pMin1, pMax1);
  //     }
  // };
  
  // template<intT dim, class nodeT>
  // struct queryBoxWard{
  //     inline tuple<point<dim>, point<dim>> getBox(nodeT* query, double r){
  //           point<dim> pMin1 = point<dim>();
  //           point<dim> pMax1 = point<dim>();
  //           for (int i=0; i<dim; ++i) {
  //           pMin1.updateX(i, query->center[i]-r);
  //           pMax1.updateX(i, query->center[i]+r);}
  //           return make_tuple(pMin1, pMax1);
  //     }
  // };
  
  // // find points in `root` within box of side `r` centered at `query` point
  // template<intT dim, class nodeT>
  // struct queryBoxAve{

  //     inline tuple<point<dim>, point<dim>> getBox(nodeT* query, double r){
  //           point<dim> pMin1 = point<dim>();
  //           point<dim> pMax1 = point<dim>();
  //           for (int i=0; i<dim; ++i) {
  //           pMin1.updateX(i, query->pMin[i]-r);
  //           pMax1.updateX(i, query->pMax[i]+r);}
  //           return make_tuple(pMin1, pMax1);
  //     }
  // };

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