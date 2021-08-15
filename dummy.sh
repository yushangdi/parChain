#!/usr/bin/env bash
echo "entering parChain/parchain/linkage/framework"
cd parchain/linkage/framework
echo "compile PC..."
make clean
GCILK=1 make -j

echo "exiting parChain/parchain/linkage/framework"
cd ../../../

./parchain/linkage/framework/linkage -method dummy -r 1 -d 2 -cachesize 32 -eps 1e-20 ./datasets/2D_GaussianDisc_10K.pbbs 