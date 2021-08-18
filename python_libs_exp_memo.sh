#!/usr/bin/env bash


[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/memo_outputs" ] || mkdir "outputs/memo_outputs"
  
methods=("ward" "avg" "complete" "avgsq")

datasets=(
			"2D_GaussianDisc_1K"
      "2D_GaussianDisc_3K"
      "2D_GaussianDisc_10K"
)
for dataset in "${datasets[@]}"; do
	for METHOD in "${methods[@]}"; do
    		echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/sklearn_${METHOD}memo_${dataset}_1th.out python3 python_benchmark_sklearn.py ${dataset} ${METHOD} > outputs/memo_outputs/sklearn_${METHOD}memo_${dataset}_1th.txt"
    		valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/sklearn_${METHOD}memo_${dataset}_1th.out python3 python_benchmark_sklearn.py ${dataset} ${METHOD} > outputs/memo_outputs/sklearn_${METHOD}memo_${dataset}_1th.txt 

				echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/scipy_${METHOD}memo_${dataset}_1th.out python3 python_benchmark_scipy.py ${dataset} ${METHOD} > outputs/memo_outputs/scipy_${METHOD}memo_${dataset}_1th.txt"
    		valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/scipy_${METHOD}memo_${dataset}_1th.out python3 python_benchmark_scipy.py ${dataset} ${METHOD} > outputs/memo_outputs/scipy_${METHOD}memo_${dataset}_1th.txt 
				
				echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/fastcluster_${METHOD}memo_${dataset}_1th.out python3 python_benchmark_fastcluster.py ${dataset} ${METHOD} > outputs/memo_outputs/fastcluster_${METHOD}memo_${dataset}_1th.txt"
				valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/fastcluster_${METHOD}memo_${dataset}_1th.out python3 python_benchmark_fastcluster.py ${dataset} ${METHOD} > outputs/memo_outputs/fastcluster_${METHOD}memo_${dataset}_1th.txt 
	
					echo "valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/fastcluster_${METHOD}-genericmemo_${dataset}_1th.out python3 python_benchmark_fastcluster.py ${dataset} ${METHOD}-generic > outputs/memo_outputs/fastcluster_${METHOD}-genericmemo_${dataset}_1th.txt"
				valgrind --tool=massif --max-snapshots=10 --detailed-freq=5 --massif-out-file=outputs/memo_outputs/fastcluster_${METHOD}-genericmemo_${dataset}_1th.out python3 python_benchmark_fastcluster.py ${dataset} ${METHOD}-generic > outputs/memo_outputs/fastcluster_${METHOD}-genericmemo_${dataset}_1th.txt 
	done
	echo "##############"
done



echo "done"