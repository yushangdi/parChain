#!/usr/bin/env bash

METHOD="$1"
SUFFIX="$2"
echo "$METHOD"


if [[ -z "$METHOD" ]]; then
  echo "Method is empty"
  exit
fi

CACHESIZE=$SUFFIX
if [[ -z "$SUFFIX" ]]; then
  SUFFIX=32
  CACHESIZE=1
fi

if [[ "avg" != "$METHOD" ]]; then
  SUFFIX=""
fi
echo "Folder suffix=$SUFFIX (only for avg-2)"  
echo "cache size=$CACHESIZE * 2 (only for avg-2)"  


[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/stats2_outputs_${METHOD}${SUFFIX}" ] || mkdir "outputs/stats2_outputs_${METHOD}${SUFFIX}"

echo "entering parChain/parchain/linkage/framework"
cd parchain/linkage/framework
echo "compile PC..."
make clean
GCILK=1 make -j

echo "exiting parChain/parchain/linkage/framework"
cd ../../../


workers=(96)

datasets=(
    "10D_UCI1_19K" 
    "2D_GaussianDisc_10K"
    "10D_UCI4_100K"
    "2D_GaussianDisc_1M"
    "5D_GaussianDisc_1M"
    "2D_UniformFill_1M"
    "5D_UniformFill_1M"
    "HT"
    "2D_GaussianDisc_10M"
    "5D_GaussianDisc_10M"
    "2D_UniformFill_10M"
    "5D_UniformFill_10M"
    "CHEM"
    "3D_GeoLife_24M"
)


        # "2D_GaussianDisc_100K"

dims=(10 2 10 2 5 2 5 10 2 5 2 5 16 3)
eps=("1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20")


for wk in "${workers[@]}"; do
    ind=0
    for dataset in "${datasets[@]}"; do
    
	if [[ "${wk}" -eq 1 ]];then
            echo "CILK_NWORKERS=${wk} ./parchain/linkage/framework/linkage -method $METHOD -r 1 -d ${dims[$ind]} -cachesize ${CACHESIZE} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs/stats2_outputs_${METHOD}${SUFFIX}/stats2_${METHOD}${SUFFIX}_${dataset}_${wk}th.txt"
            CILK_NWORKERS=${wk} ./parchain/linkage/framework/linkage -method $METHOD -r 1 -d ${dims[$ind]} -cachesize ${CACHESIZE} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                > outputs/stats2_outputs_${METHOD}${SUFFIX}/stats2_${METHOD}${SUFFIX}_${dataset}_${wk}th.txt
        else
            echo "CILK_NWORKERS=${wk} numactl -i all ./parchain/linkage/framework/linkage -method $METHOD -r 1 -cachesize ${CACHESIZE} -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs/stats2_outputs_${METHOD}${SUFFIX}/stats2_${METHOD}${SUFFIX}_${dataset}_${wk}th.txt"
            CILK_NWORKERS=${wk} numactl -i all ./parchain/linkage/framework/linkage -method $METHOD -r 1 -cachesize ${CACHESIZE} -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                > outputs/stats2_outputs_${METHOD}${SUFFIX}/stats2_${METHOD}${SUFFIX}_${dataset}_${wk}th.txt
        fi
        let ind++
    done
done

echo "done"
