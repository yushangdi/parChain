#!/usr/bin/env bash

[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/memo_outputs" ] || mkdir "outputs/memo_outputs"


echo "entering /parChain/fast_protein_cluster"
cd fast_protein_cluster
echo "compile Jeon with debug flag..."
make clean
make debug

echo "exiting /parChain/fast_protein_cluster"
cd ../

datasets=(
      "2D_GaussianDisc_1K"
      "2D_GaussianDisc_3K"
      "2D_GaussianDisc_10K"
)
for dataset in "${datasets[@]}"; do
  echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/fastprotein_avgmemo_${dataset}_36th.out fast_protein_cluster/fast_protein_cluster -i ./datasets/ ${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads 32 --haverage --ndim 2 &> outputs/memo_outputs/fastprotein_avgmemo_${dataset}_36th.txt"
  valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/fastprotein_avgmemo_${dataset}_36th.out fast_protein_cluster/fast_protein_cluster -i ./datasets/${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads 32 --haverage --ndim 2 &> outputs/memo_outputs/fastprotein_avgmemo_${dataset}_36th.txt
done

for dataset in "${datasets[@]}"; do
  echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/fastprotein_completememo_${dataset}_36th.out fast_protein_cluster/fast_protein_cluster -i ./datasets/${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads 32 --hcomplete --ndim 2 &> outputs/memo_outputs/fastprotein_completememo_${dataset}_36th.txt"
  valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/fastprotein_completememo_${dataset}_36th.out fast_protein_cluster/fast_protein_cluster -i ./datasets/${dataset}.pbbs -o output --euclidean --nclusters 1 --nthreads 32 --hcomplete --ndim 2 &> outputs/memo_outputs/fastprotein_completememo_${dataset}_36th.txt
done
echo "done"