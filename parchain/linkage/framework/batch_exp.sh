#!/usr/bin/env bash

METHOD="$1"
SUFFIX="$2"
echo "$METHOD"
echo "$SUFFIX"  # cachesize
CACHESIZE=$SUFFIX
if [[ -z "$METHOD" ]]; then
  echo "Method is empty"
  exit
fi

if [[ -z "$SUFFIX" ]]; then
  SUFFIX=32
  CACHESIZE=1
fi
    if [[ "avg" != "$METHOD" ]]; then
      SUFFIX=""
    fi
    echo "$SUFFIX" 

[ -d "outputs_${METHOD}${SUFFIX}" ] || mkdir "outputs_${METHOD}${SUFFIX}"

#make clean
#GCILK=1 make -j

workers=(96 48 36 24 12 1) #  96 72  12  48 

# datasets=(
#     "2D_GaussianDisc_100K"
#     "2D_GaussianDisc_1M"
#     "5D_GaussianDisc_1M"
#     "2D_GaussianDisc_10M"
#     "5D_GaussianDisc_10M"
# )
# #                
# dims=(2 2 2 5 2 5)



datasets=(
    "5D_UniformFill_10M"
)
    # "2D_GaussianDisc_10M"
    # "5D_GaussianDisc_10M"
    # "CHEM"
    # "3D_GeoLife_24M"
    # "5D_UniformFill_1M"
    # "2D_UniformFill_10M"
#         "5D_UniformFill_10M"         
dims=(5 10 10 2 2 2 5 10 2 5 16 3 5 2)
eps=("1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20")


for wk in "${workers[@]}"; do
    ind=0
    for dataset in "${datasets[@]}"; do
    
	if [[ "${wk}" -eq 1 ]];then
            echo "CILK_NWORKERS=${wk} ./linkage -method $METHOD -r 1 -d ${dims[$ind]} -cachesize ${CACHESIZE} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs_${METHOD}${SUFFIX}/${METHOD}${SUFFIX}_${dataset}_${wk}th.txt"
            CILK_NWORKERS=${wk} ./linkage -method $METHOD -r 1 -d ${dims[$ind]} -cachesize ${CACHESIZE} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                > outputs_${METHOD}${SUFFIX}/${METHOD}${SUFFIX}_${dataset}_${wk}th.txt
        else
            echo "CILK_NWORKERS=${wk} numactl -i all ./linkage -method $METHOD -r 3 -cachesize ${CACHESIZE} -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs_${METHOD}${SUFFIX}/${METHOD}${SUFFIX}_${dataset}_${wk}th.txt"
            CILK_NWORKERS=${wk} numactl -i all ./linkage -method $METHOD -r 3 -cachesize ${CACHESIZE} -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                > outputs_${METHOD}${SUFFIX}/${METHOD}${SUFFIX}_${dataset}_${wk}th.txt
        fi
        let ind++
    done
done

echo "done"
