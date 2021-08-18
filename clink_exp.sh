#!/usr/bin/env bash

[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/outputs_clink" ] || mkdir "outputs/outputs_clink"

# git submodule init

echo "entering parChain/althuas"
cd althuas
echo "compile Althus..."
make clean
make -j

echo "exiting parChain/althaus"
cd ../


datasets=(
    "10D_UCI1_19K" 
    "2D_GaussianDisc_10K"
    "10D_UCI4_100K"
)      
#    
dims=(10 2 10)
sizes=(19020 10000 100000)
workers=(96 48 36 24 12 1)

make clean; make;

for wk in "${workers[@]}"; do
  ind=0
  for dataset in "${datasets[@]}"; do
  	if [[ "${wk}" -eq 1 ]];then
    echo "CILK_NWORKERS=${wk} althuas/clink datasets/${dataset}.pbbs ${sizes[$ind]} ${dims[$ind]} 128 >> outputs/outputs_clink/clink_${wk}th.txt"
    echo ${dataset} >> outputs/outputs_clink/clink_${wk}th.txt
    CILK_NWORKERS=${wk} althuas/clink datasets/${dataset}.pbbs ${sizes[$ind]} ${dims[$ind]} 128 >> outputs/outputs_clink/clink_${wk}th.txt
    else
    echo "CILK_NWORKERS=${wk} numactl -i all  althuas/clink datasets/${dataset}.pbbs ${sizes[$ind]} ${dims[$ind]} 128 >> outputs/outputs_clink/clink_${wk}th.txt"
    echo ${dataset} >> outputs/outputs_clink/clink_${wk}th.txt
    CILK_NWORKERS=${wk} numactl -i all  althuas/clink datasets/${dataset}.pbbs ${sizes[$ind]} ${dims[$ind]} 128 >> outputs/outputs_clink/clink_${wk}th.txt
    fi
    let ind++
  done
done

echo "done"

[ -d "outputs/memo_outputs" ] || mkdir "outputs/memo_outputs"

echo "entering parChain/althuas"
cd althuas
echo "compile Althus..."
make clean
make debug

echo "exiting parChain/althaus"
cd ../

echo "valgrind --trace-children=yes --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/clinkmemo_2D_GaussianDisc_1K_36th.out env CILK_NWORKERS=36 althuas/clink datasets/2D_GaussianDisc_1K.pbbs 1000 2 128 >> outputs/memo_outputs/clinkmemo_2D_GaussianDisc_1K_36th.txt; "
valgrind --trace-children=yes --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/clinkmemo_2D_GaussianDisc_1K_36th.out env CILK_NWORKERS=36 althuas/clink datasets/2D_GaussianDisc_1K.pbbs 1000 2 128 >> outputs/memo_outputs/clinkmemo_2D_GaussianDisc_1K_36th.txt; 

echo "valgrind --trace-children=yes  --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/clinkmemo_2D_GaussianDisc_3K_36th.out env CILK_NWORKERS=36 althuas/clink datasets/2D_GaussianDisc_3K.pbbs 3000 2 128 >> outputs/memo_outputs/clinkmemo_2D_GaussianDisc_3K_36th.txt;"
valgrind --trace-children=yes  --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/clinkmemo_2D_GaussianDisc_3K_36th.out env CILK_NWORKERS=36 althuas/clink datasets/2D_GaussianDisc_3K.pbbs 3000 2 128 >> outputs/memo_outputs/clinkmemo_2D_GaussianDisc_3K_36th.txt;

echo "valgrind --trace-children=yes  --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/clinkmemo_2D_GaussianDisc_10K_36th.out env CILK_NWORKERS=36 althuas/clink datasets/2D_GaussianDisc_10K.pbbs 10000 2 128 >> outputs/memo_outputs/clinkmemo_2D_GaussianDisc_10K_36th.txt"
valgrind --trace-children=yes  --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/clinkmemo_2D_GaussianDisc_10K_36th.out env CILK_NWORKERS=36 althuas/clink datasets/2D_GaussianDisc_10K.pbbs 10000 2 128 >> outputs/memo_outputs/clinkmemo_2D_GaussianDisc_10K_36th.txt