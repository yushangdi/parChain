// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <algorithm>
#include "gettime.h"
#include "utils.h"
#include "geometry.h"
#include "geometryIO.h"
#include "parseCommandLine.h"
#include "parallel.h"
#include "linkage.h"
using namespace std;
using namespace benchIO;

// *************************************************************
//  TIMING
// *************************************************************

template<int dim>
void timeLinkage(point<dim>* pts, intT n, int rounds, char* outFile, commandLine P) {
  UnionFind::ParUF<int> * R;
  for (int i=0; i < rounds; i++) {
    startTime();
    UnionFind::ParUF<intT> *R = new UnionFind::ParUF<intT>(n, true);
    R = linkage<dim>(pts, n, P, R);
    R->del();
    delete R;
    nextTimeN();
  }
  cout << endl;

  // if (outFile != NULL) write(R); // todo
}

  template <class pointT>
  _seq<pointT> readPointsFromFileNoHeader(char* fname) {
    _seq<char> S = readStringFromFile(fname);
    words W = stringToWords(S.A, S.n);
    int d = pointT::dim;
    long n = (W.m)/d;
    pointT *P = newA(pointT, n);
    parsePoints(W.Strings, P, n);
    W.del();
    return _seq<pointT>(P, n);
  }

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] [-d <dim>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);

  int dim = P.getOptionIntValue("-d",2);//readPointsDimensionFromFile(iFile);

  if (dim == 1) {
    cout << "dimension 1 not supported, abort" << endl; abort();
  } else if (dim == 2) {
    _seq<point<2>> PIn = readPointsFromFileNoHeader<point<2>>(iFile);
    timeLinkage<2>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else if (dim == 3) {
    _seq<point<3>> PIn = readPointsFromFileNoHeader<point<3>>(iFile);
    timeLinkage<3>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else if (dim == 4) {
    _seq<point<4>> PIn = readPointsFromFileNoHeader<point<4>>(iFile);
    timeLinkage<4>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else if (dim == 5) {
    _seq<point<5>> PIn = readPointsFromFileNoHeader<point<5>>(iFile);
    timeLinkage<5>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else if (dim == 6) {
    _seq<point<6>> PIn = readPointsFromFileNoHeader<point<6>>(iFile);
    timeLinkage<6>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else if (dim == 7) {
    _seq<point<7>> PIn = readPointsFromFileNoHeader<point<7>>(iFile);
    timeLinkage<7>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else if (dim == 8) {
    _seq<point<8>> PIn = readPointsFromFileNoHeader<point<8>>(iFile);
    timeLinkage<8>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else if (dim == 9) {
    _seq<point<9>> PIn = readPointsFromFileNoHeader<point<9>>(iFile);
    timeLinkage<9>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  }else if (dim == 10) {
    _seq<point<10>> PIn = readPointsFromFileNoHeader<point<10>>(iFile);
    timeLinkage<10>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  }else if (dim == 16) {
    _seq<point<16>> PIn = readPointsFromFileNoHeader<point<16>>(iFile);
    timeLinkage<16>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else if (dim == 128) {
    _seq<point<128>> PIn = readPointsFromFileNoHeader<point<128>>(iFile);
    timeLinkage<128>(PIn.A, PIn.n, rounds, oFile, P);
    PIn.del();
  } else {
    cout << "dimension " << dim << " not yet supported, abort" << endl; abort();
  }
}
