#!/usr/bin/env bash

echo "1"
CILK_NWORKERS=12 numactl -i all ./linkage -method complete -r 3 -cachesize 1 -d 2 -eps 1e-20 ./datasets/2D_GaussianDisc_10M.pbbs > outputs_complete/complete_2D_GaussianDisc_10M_12th.txt

echo "2"
CILK_NWORKERS=12 numactl -i all ./linkage -method complete -r 3 -cachesize 1 -d 5 -eps 1e-20 ./datasets/5D_GaussianDisc_10M.pbbs > outputs_complete/complete_5D_GaussianDisc_10M_12th.txt

echo "3"
CILK_NWORKERS=12 numactl -i all ./linkage -method complete -r 3 -cachesize 1 -d 16 -eps 1e-20 ./datasets/CHEM.pbbs > outputs_complete/complete_CHEM_12th.txt

echo "4"
CILK_NWORKERS=24 numactl -i all ./linkage -method complete -r 3 -cachesize 1 -d 16 -eps 1e-20 ./datasets/CHEM.pbbs > outputs_complete/complete_CHEM_24th.txt

echo "5"
CILK_NWORKERS=36 numactl -i all ./linkage -method complete -r 3 -cachesize 1 -d 16 -eps 1e-20 ./datasets/CHEM.pbbs > outputs_complete/complete_CHEM_36th.txt
echo "done"
