#!/usr/bin/env bash

[ -d "outputs_matrix" ] || mkdir "outputs_matrix"

# make clean
# GCILK=1 make -j

workers=(24) # 96 48 36 24 12 1


datasets=(
    "10D_UCI4_100K" 
    "10D_UCI1_19K" 
    "2D_GaussianDisc_10K"
)

#         "5D_UniformFill_10M"         
dims=(10 10 2)
eps=("1e-20" "1e-20" "1e-20" )

methods=("ward" "avg" "complete" "avgsq") # 
for METHOD in "${methods[@]}"; do
    for wk in "${workers[@]}"; do
        ind=0
        for dataset in "${datasets[@]}"; do
            if [[ "${wk}" -eq 1 ]];then
                echo "CILK_NWORKERS=${wk} ./linkage -method $METHOD -r 1 -d ${dims[$ind]} -matrix -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs_matrix/${METHOD}_matrix_${dataset}_${wk}th.txt"
                CILK_NWORKERS=${wk} ./linkage -method $METHOD -r 1 -d ${dims[$ind]} -matrix -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                    > outputs_matrix/${METHOD}_matrix_${dataset}_${wk}th.txt
            else
                echo "CILK_NWORKERS=${wk} numactl -i all ./linkage -method $METHOD -r 3 -matrix -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > outputs_matrix/${METHOD}_matrix_${dataset}_${wk}th.txt"
                CILK_NWORKERS=${wk} numactl -i all ./linkage -method $METHOD -r 3 -matrix -d ${dims[$ind]} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                    > outputs_matrix/${METHOD}_matrix_${dataset}_${wk}th.txt
            fi
            let ind++
        done
    done
done

echo "done"