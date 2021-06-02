#!/usr/bin/env bash

METHOD="$1"
echo "$METHOD"

if [[ -z "$METHOD" ]]; then
  echo "Method is empty"
  exit
fi

[ -d "outputs_opt" ] || mkdir "outputs_opt"

# make clean
# GCILK=1 make -j

workers=(96) #     

datasets=(
    "10D_UCI4_100K" 
    "2D_UniformFill_1M"
)
dims=(10 2)
eps=("1e-20" "1e-20")

for wk in "${workers[@]}"; do
    ind=0
    for dataset in "${datasets[@]}"; do
        if [[ "${wk}" -eq 1 ]];then
            echo "CILK_NWORKERS=${wk} ./linkage -naivethresh 100000000 -method $METHOD -r 1 -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs_${METHOD}/${METHOD}_${dataset}_${wk}th.txt"
            CILK_NWORKERS=${wk} ./linkage -naivethresh  100000000 -method $METHOD -r 1 -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                > outputs_opt/${METHOD}_norange_${dataset}_${wk}th.txt
        else
            echo "CILK_NWORKERS=${wk} numactl -i all ./linkage  -naivethresh  100000000 -method $METHOD -r 3 -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs_${METHOD}/${METHOD}_${dataset}_${wk}th.txt"
            CILK_NWORKERS=${wk} numactl -i all ./linkage  -naivethresh  100000000 -method $METHOD -r 1 -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                > outputs_opt/${METHOD}_norange_${dataset}_${wk}th.txt
        fi
        let ind++
    done
done

echo "done"