# ParChain: A Framework for Parallel Hierarchical Agglomerative Clustering using Nearest-Neighbor Chain

This program is adapted from [the problem based benchmark suite (PBBS)](https://github.mit.edu/jshun/pbbs).


Implementations
--------

`parchain`
Our ParChain algorithms for complete linkage with Euclidean distance, Ward's linkage with Euclidean distance, average linkage with Euclidean distance and average linkage with squared Euclidean distance .

`multi-threaded-NN-chain`
[Jeon et al's](https://ieeexplore.ieee.org/document/6893001) parallel nearest-neighbor chain implementation. It supports average linkage with Euclidean distance. 

`fast_protein_cluster`
[fast_protein_cluster](https://pubmed.ncbi.nlm.nih.gov/24532722/) parallel HAC implementation that works for complete and average linkage with Euclidean distance metric.
The [original code](https://github.com/lhhunghimself/fast_protein_cluster) works for root-mean-square deviation (RMSD) and template modeling score (TM-score) measures. 


## Installation

Compiler:
* g++ = 7.5.0 with support for Cilk Plus

***********************************
TEST DATA
***********************************

/data

