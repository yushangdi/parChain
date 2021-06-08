#!/usr/bin/env bash


[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/memo_outputs" ] || mkdir "outputs/memo_outputs"

echo "entering parChain/parchain/linkage/framework"
cd parchain/linkage/framework
echo "compile PC..."
make clean
GCILKDEBUG=1 make -j

echo "exiting parChain/parchain/linkage/framework"
cd ../../../

methods=("ward" "avg" "complete" "avgsq")

datasets=(
    "2D_GaussianDisc_1K"
    "2D_GaussianDisc_3K"
    "2D_GaussianDisc_10K"
)

dims=(2 2 2)
eps=("1e-20" "1e-20" "1e-20")


ind=0
for dataset in "${datasets[@]}"; do
  for METHOD in "${methods[@]}"; do
    SUFFIX=32
    CACHESIZE=32
    if [[ "avg" != "$METHOD" ]]; then
      SUFFIX=""
      CACHESIZE=1
    fi
    echo "$SUFFIX" 
    echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.out ./parchain/linkage/framework/linkage -method $METHOD -r 1 -d ${dims[$ind]} -cachesize ${CACHESIZE} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs/memo_outputs/${METHOD}${SUFFIX}_${dataset}_96th.txt"
    valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.out ./parchain/linkage/framework/linkage -method $METHOD -r 1 -d ${dims[$ind]} -cachesize ${CACHESIZE} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
            > outputs/memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.txt
  done
  let ind++
done

echo "done"
