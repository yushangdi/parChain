CC = g++
CCFLAGS = -O3 -Wall -std=c++17

NNOBJECTS =  cscipy.o

cscipy: $(NNOBJECTS) Makefile
	$(CC) $(CCFLAGS) $(NNOBJECTS) -o $@

cscipy.o: cscipy.cpp
	$(CC) $(CCFLAGS) -c $<

debug: CCFLAGS= -g -Wall -std=c++17
 
debug: cscipy

clean:
	rm -f cscipy *.o
