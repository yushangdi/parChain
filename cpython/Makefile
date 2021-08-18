CC = g++
CCFLAGS = -O3 -Wall -std=c++17

NNOBJECTS =  cpython.o

cpython: $(NNOBJECTS) Makefile
	$(CC) $(CCFLAGS) $(NNOBJECTS) -o $@

cpython.o: cpython.cpp
	$(CC) $(CCFLAGS) -c $<

debug: CCFLAGS= -g -Wall -std=c++17
 
debug: cpython

clean:
	rm -f cpython *.o
