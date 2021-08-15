#!/usr/bin/env bash
echo "entering parChain/parchain/linkage/framework"
cd parchain/linkage/framework
echo "compile PC..."
make clean
GCILK=1 make -j

echo "exiting parChain/parchain/linkage/framework"
cd ../../../

METHOD="dummy"
SUFFIX="32"

[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/outputs_${METHOD}${SUFFIX}" ] || mkdir "outputs/outputs_${METHOD}${SUFFIX}"

numactl -i all ./parchain/linkage/framework/linkage -method dummy -r 1 -d 2 -cachesize 32 -eps 1e-20 ./datasets/2D_GaussianDisc_10K.pbbs  > outputs/outputs_${METHOD}${SUFFIX}/${METHOD}${SUFFIX}_2D_GaussianDisc_10K_96th.txt
numactl -i all ./parchain/linkage/framework/linkage -method dummy -r 1 -d 10 -cachesize 32 -eps 1e-20 ./datasets/10D_UCI1_19K.pbbs  > outputs/outputs_${METHOD}${SUFFIX}/${METHOD}${SUFFIX}_10D_UCI1_19K_96th.txt
numactl -i all ./parchain/linkage/framework/linkage -method dummy -r 1 -d 10 -cachesize 32 -eps 1e-20 ./datasets/10D_UCI4_100K.pbbs  > outputs/outputs_${METHOD}${SUFFIX}/${METHOD}${SUFFIX}_10D_UCI4_100K_96th.txt
