include parallelDefs

# required files
NEAREST_NEIGHBORS = geometry.h transpose.h
GLOBAL_REQUIRE = parallel.h sequence.h $(NEAREST_NEIGHBORS) sampleSort.h
LOCAL_REQUIRE = unionfind.h
OBJS = framework.o

include ../common/timeRequiredFiles
include MakeBench
