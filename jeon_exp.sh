#!/usr/bin/env bash

[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/outputs_jeon" ] || mkdir "outputs/outputs_jeon"

echo "entering parChain/multi-threaded-NN-chain"
cd multi-threaded-NN-chain
echo "compile Jeon..."
make clean
make -j

echo "exiting parChain/multi-threaded-NN-chain"
cd ../

workers=(96 48 36 24 12 8 4 2 1) 
workers2=(95 47 35 23 11 7 3 1 1)  

# workers=(1)
# workers2=(1)

# workers used for distance matrix computations, chosen for each dataset 
# to optimize speed
workers3=(12 12 24)
datasets=(
    "10D_UCI1_19K" 
    "2D_GaussianDisc_10K"
    "10D_UCI4_100K"
)      
#    
dims=(10 2 10)
sizes=(19020 10000 100000)

j=0
for wk in "${workers2[@]}"; do
    ind=0
    numa="numactl -C 1 "
    wk3=${workers3[$ind]}
    for dataset in "${datasets[@]}"; do
        if [[ "${workers[$j]}" -eq 1 ]];then
        	rd=1
        else
        	rd=3
        	numa="numactl -i all "
	fi
	if [[ "${wk3}" -ge ${workers[$j]} ]];then
                wk3=${workers[$j]}
        fi

        echo "${numa}./multi-threaded-NN-chain/NN-chain ./datasets/${dataset}.pbbs ${sizes[$ind]} ${dims[$ind]} ${wk} 1 ${wk3} ${rd} > outputs/outputs_jeon/jeon_${dataset}_${workers[$j]}th.txt"
        ${numa}./multi-threaded-NN-chain/NN-chain ./datasets/${dataset}.pbbs ${sizes[$ind]} ${dims[$ind]} ${wk} 1 ${wk3} ${rd} > outputs/outputs_jeon/jeon_${dataset}_${workers[$j]}th.txt
        let ind++
	echo
    done
    echo "######"
    let j++
done

echo "done"