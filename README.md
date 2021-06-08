# ParChain: A Framework for Parallel Hierarchical Agglomerative Clustering using Nearest-Neighbor Chain

This program is adapted from [the problem based benchmark suite (PBBS)](https://github.mit.edu/jshun/pbbs).

To get the submodules:
```bash
git pull
git submodule update --init
```

Implementations
--------
Each folder contains an implementation that we tested.

`parchain`
Our ParChain algorithms for complete linkage with Euclidean distance, Ward's linkage with Euclidean distance, average linkage with Euclidean distance and average linkage with squared Euclidean distance .

`multi-threaded-NN-chain`
[Jeon et al's](https://ieeexplore.ieee.org/document/6893001) parallel nearest-neighbor chain implementation. It supports average linkage with Euclidean distance. 

`fast_protein_cluster`
[fast_protein_cluster](https://pubmed.ncbi.nlm.nih.gov/24532722/) parallel HAC implementation that works for complete and average linkage with Euclidean distance metric.
The [original code](https://github.com/lhhunghimself/fast_protein_cluster) works for root-mean-square deviation (RMSD) and template modeling score (TM-score) measures. 

`fastcluster`
[fastcluster](http://danifold.net/fastcluster.html?section=1) is a sequential HAC implementation in C++ with python interface. 

`scipy` and `sklearn` can be installed from python.

## Installation

Compiler:
* g++ = 7.5.0 with support for Cilk Plus
* python3 &gt;= 3.6.9
* make

Python libraries:
* scipy &gt;= 1.5.4 
* sklearn &gt;= 0.24.2
* pandas
* numpy

## Data

You can download our datasets [here](https://console.cloud.google.com/storage/browser/...)


## Running Tests
For running time tests, we use `numactl`. It can be installed using `apt install numactl`. 
For the memory tests, `valgrind &gt;= 3.17.0` is required. You can download [here](https://www.valgrind.org/docs/download_docs.html). 

To run our scripts, the datasets need to be in the "datasets" folder under this folder. 

### PC and PC-matrix

The runtime/memory experiment result of our algorithm runs the following. The result files will be in `outputs/outputs_{method}` and `outputs/memo_outputs`.
You can modify `datasets` and `dims` in pc_exp.sh for the datasets to run. 
```bash
./pc_exp.sh complete
./pc_exp.sh ward
./pc_exp.sh avg 32
./pc_exp.sh avgsq 

./pc_exp_memo.sh
```

The runtime/memory experiment result of our algorithm runs the following. The result files will be in `outputs/outputs_matrix` and `outputs/memo_outputs`.
You can modify `datasets` and `dims` in pc_exp.sh for the datasets to run. 
```bash
./pc_matrix_exp.sh complete
./pc_matrix_exp.sh ward
./pc_matrix_exp.sh avg
./pc_matrix_exp.sh avgsq 

./pc_matrix_exp_memo.sh
```
 

### Jeon
The runtime/memory experiment result of Jeon runs the following. The result files will be in `outputs/outputs_jeon` and `outputs/memo_outputs`.
```bash
./jeon_exp.sh
./jeon_exp_memo.sh
```

### fastprotein
The runtime/memory experiment result of fastprotein runs the following. The result files will be in `outputs/outputs_fastprotein` and `outputs/memo_outputs`.
```bash
./fastprotein_exp.sh
./fastprotein_exp_memo.sh
```


### Scipy and Sklearn and fastcluster
The runtime/memory experiment result of scipy and sklearn runs the following. The result files will be in `outputs/outputs_scipy`, `outputs/outputs_sklearn`, `outputs/outputs_fastcluster` and `outputs/memo_outputs`.
You need to install Scipy, Sklearn, and fastcluster first. 
Scipy and Sklearn and be installed using `pip`. The fastcluster submodule in this repo includes some additional methods (linear memory nearest neighbor chain for Ward's and average linkage with squared Euclidean distance metric). You can refer to `fastcluster/` for installation method. 
```bash
./python_libs_exp.sh
./python_libs_exp_memo.sh
```