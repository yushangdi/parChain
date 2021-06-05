#!/usr/bin/env bash

[ -d "memo_outputs" ] || mkdir "memo_outputs"
  
# methods=("avg" "complete")
#methods=("avg")
datasets=(
      "2D_GaussianDisc_10K"
      "2D_GaussianDisc_100K"
)
for dataset in "${datasets[@]}"; do
  echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=memo_outputs/fastprotein_avgmemo_${dataset}_36th.out ./fast_protein_cluster -i ~/datasets/ ${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads 32 --haverage --ndim 2 > memo_outputs/fastprotein_avgmemo_${dataset}_36th.txt"
  valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=memo_outputs/fastprotein_avgmemo_${dataset}_36th.out ./fast_protein_cluster -i ~/datasets/ ${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads 32 --haverage --ndim 2 > memo_outputs/fastprotein_avgmemo_${dataset}_36th.txt
done

for dataset in "${datasets[@]}"; do
  echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=memo_outputs/fastprotein_completememo_${dataset}_36th.out ./fast_protein_cluster -i ~/datasets/ ${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads 32 --hcomplete --ndim 2 > memo_outputs/fastprotein_completememo_${dataset}_36th.txt"
  valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=memo_outputs/fastprotein_completememo_${dataset}_36th.out ./fast_protein_cluster -i ~/datasets/ ${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads 32 --hcomplete --ndim 2 > memo_outputs/fastprotein_completememo_${dataset}_36th.txt
done
echo "done"