#!/usr/bin/env bash


[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/memo_outputs" ] || mkdir "outputs/memo_outputs"

echo "entering ~/parChain/multi-threaded-NN-chain"
cd ~/parChain/multi-threaded-NN-chain
echo "compile Jeon with debug flag..."
make clean
make debug

echo "exiting ~/parChain/multi-threaded-NN-chain"
cd ../

echo "valgrind --tool=massif --massif-out-file=outputs/memo_outputs/jeonmemo_2D_GaussianDisc_1K_4th.out ./multi-threaded-NN-chain/NN-chain ./datasets/2D_GaussianDisc_1K.pbbs 1000 2 3 1 4 1 > outputs/memo_outputs/jeonmemo_2D_GaussianDisc_1K_4th.txt"
valgrind --tool=massif --massif-out-file=outputs/memo_outputs/jeonmemo_2D_GaussianDisc_1K_4th.out ./multi-threaded-NN-chain/NN-chain ./datasets/2D_GaussianDisc_1K.pbbs 1000 2 3 1 4 1 > outputs/memo_outputs/jeonmemo_2D_GaussianDisc_1K_4th.txt

echo "valgrind --tool=massif --massif-out-file=outputs/memo_outputs/jeonmemo_2D_GaussianDisc_5K_4th.out ./multi-threaded-NN-chain/NN-chain ./datasets/2D_GaussianDisc_5K.pbbs 5000 2 3 1 4 1 > outputs/memo_outputs/jeonmemo_2D_GaussianDisc_5K_4th.txt"
valgrind --tool=massif --massif-out-file=outputs/memo_outputs/jeonmemo_2D_GaussianDisc_5K_4th.out ./multi-threaded-NN-chain/NN-chain ./datasets/2D_GaussianDisc_5K.pbbs 5000 2 3 1 4 1 > outputs/memo_outputs/jeonmemo_2D_GaussianDisc_5K_4th.txt

echo "valgrind --tool=massif --massif-out-file=outputs/memo_outputs/jeonmemo_2D_GaussianDisc_100K_4th.out ./multi-threaded-NN-chain/NN-chain ./datasets/2D_GaussianDisc_10K.pbbs 10000 2 3 1 4 1 > outputs/memo_outputs/jeonmemo_2D_GaussianDisc_10K_4th.txt"
valgrind --tool=massif --massif-out-file=outputs/memo_outputs/jeonmemo_2D_GaussianDisc_100K_4th.out ./multi-threaded-NN-chain/NN-chain ./datasets/2D_GaussianDisc_10K.pbbs 10000 2 3 1 4 1 > outputs/memo_outputs/jeonmemo_2D_GaussianDisc_10K_4th.txt

echo "done"
