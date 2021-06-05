#!/usr/bin/env bash


CILK_NWORKERS=96 numactl -i all ./linkage -method avg -r 3 -d 2 -eps 1e-20 -cachesize 256 ./datasets/2D_GaussianDisc_10M.pbbs > debug/avg256_2D_GaussianDisc_10M_96th.txt
CILK_NWORKERS=96 numactl -i all ./linkage -method avg -r 3 -d 2 -eps 1e-20 -cachesize 128 ./datasets/2D_GaussianDisc_10M.pbbs > debug/avg_2D_GaussianDisc_10M_96th.txt
CILK_NWORKERS=96 numactl -i all ./linkage -method avg -r 3 -d 2 -eps 1e-20 -cachesize 64 ./datasets/2D_GaussianDisc_10M.pbbs > debnug/avg64_2D_GaussianDisc_10M_96th.txt


# python3 python_benchmark.py 10D_UCI1_19K all
# python3 python_benchmark.py 10D_UCI4_100K all
# python3 python_benchmark.py 2D_GaussianDisc_10K all

# python3 python_benchmark_sklearn.py 10D_UCI1_19K all
# python3 python_benchmark_sklearn.py 10D_UCI4_100K all
# python3 python_benchmark_sklearn.py 2D_GaussianDisc_10K all


# #  sudo python3 setup.py install
# python3 python_benchmark_fastcluster.py 10D_UCI1_19K all 0
# python3 python_benchmark_fastcluster.py 10D_UCI4_100K all 0
# python3 python_benchmark_fastcluster.py 2D_GaussianDisc_10K all 0
# python3 python_benchmark_fastcluster.py 2D_UniformFill_1M linear 1
# python3 python_benchmark_fastcluster.py 2D_VisualVar_1M linear 1
# python3 python_benchmark_fastcluster.py 5D_UniformFill_1M linear 1
# python3 python_benchmark_fastcluster.py 5D_VisualVar_1M linear 1
# python3 python_benchmark_fastcluster.py HT linear 1

# echo "python3 python_benchmark_fastcluster.py 2D_VisualVar_1M linear, 2hr limit"
# python3 python_benchmark_fastcluster.py 2D_VisualVar_1M linear 1 &
# pid=$!
# sleep 7200 # two hours
# if kill -0 $pid; then kill -TERM $pid; echo "too long"; fi


# jeon et al
# change pdist nnumthreads to 1 first
# can't change numthreds in nnchain to 1, will not terminate
# numactl -C 1 ./NN-chain ~/datasets/2D_GaussianDisc_10K.pbbs 10000 2 1 1 1 3
# ./NN-chain ~/datasets/2D_GaussianDisc_10K.pbbs 10000 2 1 1 2 3
# ./NN-chain ~/datasets/2D_GaussianDisc_10K.pbbs 10000 2 3 1 4 3
# ./NN-chain ~/datasets/2D_GaussianDisc_10K.pbbs 10000 2 7 1 8 3
# ./NN-chain ~/datasets/2D_GaussianDisc_10K.pbbs 10000 2 15 1 16 3
# ./NN-chain ~/datasets/2D_GaussianDisc_10K.pbbs 10000 2 31 1 32 3


# numactl -C 1 ./NN-chain ~/datasets/10D_UCI1_19K.pbbs 19020 10 1 1 3
# ./NN-chain ~/datasets/10D_UCI1_19K.pbbs 19020 10 1 1 2 3
# ./NN-chain ~/datasets/10D_UCI1_19K.pbbs 19020 10 3 1 4 3
# ./NN-chain ~/datasets/10D_UCI1_19K.pbbs 19020 10 7 1 8 3
# ./NN-chain ~/datasets/10D_UCI1_19K.pbbs 19020 10 15 1 16 3
# ./NN-chain ~/datasets/10D_UCI1_19K.pbbs 19020 10 31 1 32 3


#proten cluster
# ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 1 --hcomplete --ndim 2
# ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 2 --hcomplete --ndim 2
# ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 4 --hcomplete --ndim 2
# ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 8 --hcomplete --ndim 2
# ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 16 --hcomplete --ndim 2
# ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 32 --hcomplete --ndim 2

# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 1 --haverage --ndim 2; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 2 --haverage --ndim 2; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 4 --haverage --ndim 2; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 8 --haverage --ndim 2; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 16 --haverage --ndim 2; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/2D_GaussianDisc_10K.pbbs -o output --euclidean --nclusters 10000 --nthreads 32 --haverage --ndim 2; echo "\n"; done

# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 1 --hcomplete --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 2 --hcomplete --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 4 --hcomplete --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 8 --hcomplete --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 16 --hcomplete --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 32 --hcomplete --ndim 10; echo "\n"; done


# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 1 --haverage --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 2 --haverage --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 4 --haverage --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 8 --haverage --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 16 --haverage --ndim 10; echo "\n"; done
# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI1_19K.pbbs -o output --euclidean --nclusters 19020 --nthreads 32 --haverage --ndim 10; echo "\n"; done


# for i in {1..3}; do  ./fast_protein_cluster -i ~/datasets/10D_UCI4_100K.pbbs -o output --euclidean --nclusters 100000 --nthreads 1 --haverage --ndim 10; echo "\n"; done
