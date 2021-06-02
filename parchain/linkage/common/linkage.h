#include "geometry.h"
#include "unionfind.h"
#include "parseCommandLine.h"

template<int dim>
UnionFind::ParUF<int> *linkage(point<dim>*, intT, commandLine, UnionFind::ParUF<intT> *);
