#!/usr/bin/env bash

[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/outputs_fastprotein" ] || mkdir "outputs/outputs_fastprotein"



echo "entering /parChain/fast_protein_cluster"
cd fast_protein_cluster
echo "compile fastprotein..."
make clean
make -j

echo "exiting /parChain/fast_protein_cluster"
cd ../

workers=(96 48 36 24 12 8 4 2 1)

datasets=(
    "10D_UCI1_19K" 
    "2D_GaussianDisc_10K"
    "10D_UCI4_100K"
)      
#    
dims=(10 2 10)

for wk in "${workers[@]}"; do
  ind=0
  for dataset in "${datasets[@]}"; do
    echo "fast_protein_cluster/fast_protein_cluster -i ./datasets/${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads ${wk} --haverage --ndim ${dims[$ind]} &> outputs/outputs_fastprotein/fastprotein_avg _${dataset}_${wk}th.txt"
    fast_protein_cluster/fast_protein_cluster -i ./datasets/${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads ${wk} --haverage --ndim ${dims[$ind]} &> outputs/outputs_fastprotein/fastprotein_avg_${dataset}_${wk}th.txt
    let ind++
  done
done

for wk in "${workers[@]}"; do
  ind=0
  for dataset in "${datasets[@]}"; do
    echo "fast_protein_cluster/fast_protein_cluster -i ./datasets/${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads ${wk} --hcomplete --ndim ${dims[$ind]} &> outputs/outputs_fastprotein/fastprotein_complete _${dataset}_${wk}th.txt"
    fast_protein_cluster/fast_protein_cluster -i ./datasets/${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads ${wk} --hcomplete --ndim ${dims[$ind]} &> outputs/outputs_fastprotein/fastprotein_complete _${dataset}_${wk}th.txt
    let ind++
  done
done 

echo "done"