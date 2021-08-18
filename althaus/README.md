# Althaus

use "git submodule init" and "git submodule update --init" to init the parlaylib submodule. 

## Compile

```bash
make 
```

or 

```bash
make debug # compile with -g flag instead of -O3
```


## Run

CILK_NWORKERS=[num. of threads] numactl -i all [dataset] [size] [dim] [num. of slots for each cluster]

example:
```bash
CILK_NWORKERS=96 numactl -i all ./datasets/2D_GaussianDisc_1K.pbbs 1000 2 128
```
