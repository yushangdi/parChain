This program is adapted from [the problem based benchmark suite (PBBS)](https://github.mit.edu/jshun/pbbs).

## Installation

Compiler:
* g++ = 7.5.0 with support for Cilk Plus

## Compile

Run the following command in the `linkage/framework` directory
```bash
GCILK=1 make
```

## Run

CILK_NWORKERS=[# of threads] numactl -i all ./linkage/framework/linkage -method [METHOD] -cachesize [half of cache size] -d [dim] [dataset]

* `METHOD` can be "complete",  "ward", "avg" (average linkage with Euclidean distance metric), or "avgsq" (average linkage with squared Euclidean distance metric).
* `cache size` is the half of the size of each hash table. if `cache size=1`, no cache will be used. It is only used for method "avg". Other methods do not keep any cache.
* Use `-r [# of times to run] ` to run the program multiple times without re-loading the dataset from the disk file.
* Use `-naivethresh [threshold]` limit the range query optimization only when the number of active clusters is larger than [threshold]. Default to 5.
* Use `-matrix` flag to use distace matrix