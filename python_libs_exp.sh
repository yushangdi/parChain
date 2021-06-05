#!/usr/bin/env bash

[ -d "outputs" ] || mkdir "outputs"
[ -d "outputs/outputs_scipy" ] || mkdir "outputs/outputs_scipy"


echo "python3 python_benchmark_scipy.py 10D_UCI1_19K all > outputs/outputs_scipy/scipy_10D_UCI1_19K_1th.txt"
python3 python_benchmark_scipy.py 10D_UCI1_19K all > outputs/outputs_scipy/scipy_10D_UCI1_19K_1th.txt

echo "python3 python_benchmark_scipy.py 10D_UCI4_100K all > outputs/outputs_scipy/scipy_10D_UCI4_100K_1th.txt"
python3 python_benchmark_scipy.py 10D_UCI4_100K all > outputs/outputs_scipy/scipy_10D_UCI4_100K_1th.txt

echo "python3 python_benchmark_scipy.py 2D_GaussianDisc_10K all > outputs/outputs_scipy/scipy_2D_GaussianDisc_10K_1th.txt"
python3 python_benchmark_scipy.py 2D_GaussianDisc_10K all > outputs/outputs_scipy/scipy_2D_GaussianDisc_10K_1th.txt

[ -d "outputs/outputs_sklearn" ] || mkdir "outputs/outputs_sklearn"

echo "python3 python_benchmark_sklearn.py 10D_UCI1_19K all > outputs/outputs_sklearn/sklearn_10D_UCI1_19K_1th.txt"
python3 python_benchmark_sklearn.py 10D_UCI1_19K all > outputs/outputs_sklearn/sklearn_10D_UCI1_19K_1th.txt

echo "python3 python_benchmark_sklearn.py 10D_UCI4_100K all > outputs/outputs_sklearn/sklearn_10D_UCI4_100K_1th.txt"
python3 python_benchmark_sklearn.py 10D_UCI4_100K all > outputs/outputs_sklearn/sklearn_10D_UCI4_100K_1th.txt

echo "python3 python_benchmark_sklearn.py 2D_GaussianDisc_10K all > outputs/outputs_sklearn/sklearn_2D_GaussianDisc_10K_1th.txt"
python3 python_benchmark_sklearn.py 2D_GaussianDisc_10K all > outputs/outputs_sklearn/sklearn_2D_GaussianDisc_10K_1th.txt

echo "done"