#!/usr/bin/env bash

[ -d "memo_outputs" ] || mkdir "memo_outputs"

# make clean
# GCILK=1 make -j

    # "2D_GaussianDisc_1K"
    # "2D_GaussianDisc_10K"

datasets=(
    "2D_GaussianDisc_3K"
)
      
dims=(2 2 2)
eps=("1e-20" "1e-20" "1e-20" )

methods=("ward" ) # "avg" "complete" "avgsq"
for METHOD in "${methods[@]}"; do
        ind=0
        for dataset in "${datasets[@]}"; do
                echo "/usr/local/bin/valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=memo_outputs/${METHOD}_matrixmemo_${dataset}_96th.out  ./linkage -method $METHOD -r 1 -d ${dims[$ind]} -matrix -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > memo_outputs/${METHOD}_matrix_${dataset}_96th.txt"
                /usr/local/bin/valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=memo_outputs/${METHOD}_matrixmemo_${dataset}_96th.out  ./linkage -method $METHOD -r 1 -d ${dims[$ind]} -matrix -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
                    > memo_outputs/${METHOD}_matrix_${dataset}_96th.txt
            let ind++
        done
done

echo "done"