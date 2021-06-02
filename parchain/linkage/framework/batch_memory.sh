#!/usr/bin/env bash


# echo "$SUFFIX"  # cachesize
CACHESIZE=32
methods=("ward" "avg" "complete" "avgsq")



# exit
[ -d "memo_outputs" ] || mkdir "memo_outputs"
    # "2D_GaussianDisc_10K"
datasets=(
    "2D_GaussianDisc_1K"
    "2D_GaussianDisc_3K"
)      
dims=(2 2)
eps=("1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20" "1e-20")


make clean
GCILK=1 make -j

ind=0
for dataset in "${datasets[@]}"; do
  for METHOD in "${methods[@]}"; do
    SUFFIX=32
    
    if [[ "avg" != "$METHOD" ]]; then
      SUFFIX=""
    fi
    echo "$SUFFIX" 
    echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.out ./linkage -method $METHOD -r 1 -d ${dims[$ind]} -cachesize ${CACHESIZE} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs > memo_outputs/${METHOD}${SUFFIX}_${dataset}_96th.txt"
    valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.out ./linkage -method $METHOD -r 1 -d ${dims[$ind]} -cachesize ${CACHESIZE} -eps ${eps[$ind]} ./datasets/${dataset}.pbbs \
            > memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.txt
  done
  let ind++
done

# echo "done"

    # echo "mem_heap_B: " >> memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.txt
    # echo "grep mem_heap_B memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 >> memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.txt"
    # grep mem_heap_B memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 >> memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.txt
    # grep mem_heap_B memo_outputs/${METHOD}${SUFFIX}memo_${dataset}_96th.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1
