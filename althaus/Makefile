CC = g++
CCFLAGS = -O3  -Wall -std=c++17 -fcilkplus -DCILKP -DPARLAY_CILK
PLFLAGS = -ldl -fcilkplus $(SDSLLF)
PARLAY_DIR=parlaylib/include
INCLUDE = -I$(PARLAY_DIR)

NNOBJECTS =  clustering.o

clink: $(NNOBJECTS) Makefile
	$(CC) $(PLFLAGS) $(NNOBJECTS) -o $@

clustering.o: clustering.cc
	$(CC) $(INCLUDE) $(CCFLAGS) -c $<

debug: CCFLAGS= -g -Wall -std=c++17  -fcilkplus -DPARLAY_USE_STD_ALLOC -DCILKP -DPARLAY_CILK

debug: clink

clean:
	rm -f clink *.o
