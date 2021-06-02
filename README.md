# ParChain: A Framework for Parallel Hierarchical Agglomerative Clustering using Nearest-Neighbor Chain

This program is adapted from [the problem based benchmark suite (PBBS)](https://github.mit.edu/jshun/pbbs).


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
* python3

Python libraries:
* scipy
* sklearn
* pandas
* numpy

***********************************
TEST DATA
***********************************

/data


## Running Tests
For the memory tests, `valgrind &gt;= 3.17.0` is required.
