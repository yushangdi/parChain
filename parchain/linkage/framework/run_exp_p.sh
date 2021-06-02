#!/bin/bash

# $1: file name
# $2: dimension of dataset
# $3: 1 if want to re compile
# $4: number of workers

if [ $3 -eq 1 ]
then
    make clean
    GCILK=1 make -j
fi

# if i > 1
echo "CILK_NWORKERS=$4 numactl -i all ./linkage -r 1 -d $2 -method $5 ./datasets/$1.pbbs"
CILK_NWORKERS=$4 numactl -i all ./linkage -r 3 -d $2 -method $5 ./datasets/$1.pbbs

# USEJEMALLOC=1  

# ./run_exp_p.sh 10D_UCI1_19K 10 1 72 complete

# 0.392
# 0.493
# 0.46...