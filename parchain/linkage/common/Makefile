include parallelDefs

ifndef USEJEMALLOC
JEMALLOC =
else
JEMALLOCLD = $(shell jemalloc-config --libdir)
JEMALLOC = -DUSEJEMALLOC -L$(JEMALLOCLD) -ljemalloc
endif

COMMON = IO.h parseCommandLine.h parallel.h runTests.py sequence.h utils.h geometry.h geometryIO.h gettime.h sampleSort.h

%.o : %.C $(COMMON)
	$(PCC) $(PCFLAGS) $(JEMALLOC) -c $< -o $@

CHECKFILES = linkageCheck.o

linkageCheck : $(CHECKFILES)
	$(PCC) $(PLFLAGS) -o $@ $(CHECKFILES)

$(COMMON) :
	ln -s ../../common/$@ .

clean :
	rm -f linkageCheck *.o

cleansrc :
	make -s clean
	rm -f $(COMMON) *.pyc
