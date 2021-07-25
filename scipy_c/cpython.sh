#!/usr/bin/env bash

[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/outputs_cscipy" ] || mkdir "outputs/outputs_cscipy"
[ -d "outputs/outputs_csklearn" ] || mkdir "outputs/outputs_csklearn"



datasets=(
    "10D_UCI1_19K" 
    "2D_GaussianDisc_10K"
    "10D_UCI4_100K"
)      
#    
dims=(10 2 10)
sizes=(19020 10000 100000)
methods={"complete", "ward", "avg", "avgsq"}
libraries={"scipy", "sklearn"}

ind=0


for dataset in "${datasets[@]}"; do
  for method in "${methods[@]}"; do
    for library in "${libraries[@]}"; do
      echo "./clink datasets/${dataset}.pbbs ${sizes[$ind]} ${dims[$ind]} 1 ${method} ${library}  >> outputs/outputs_c${library}/c${library}_${dataset}_1th.txt"
      echo ${dataset}_${method} >> outputs/outputs_c${library}/c${library}_${dataset}_1th.txt
      ./clink datasets/${dataset}.pbbs ${sizes[$ind]} ${dims[$ind]} 1 ${method} ${library}>> outputs/outputs_c${library}/c${library}_${dataset}_1th.txt
    done
  done
  let ind++
done

echo "done"