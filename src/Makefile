
CXX      = g++ -std=c++11 #-g #-pg #-fpermissive #-g # -pg
LIBGAB   = ../lib/libgab/

CXXFLAGS = -Wall -lm -O3 -lz -I../lib/libgab/ -I../lib/libgab//gzstream/ -I../lib/samtools/ -I../lib/htslib/    -c
LDFLAGS  = -lgsl -lcurl -lgslcblas -llzma -lbz2 -lz  -lpthread 
LDLIBS   = ../lib/samtools/bedidx.o ../lib/htslib/libhts.a ../lib/samtools/libbam.a ../lib/samtools/libst.a 


all:  ../lib/libgab/utils.o  ../lib/htslib/libhts.a ../lib/samtools/bedidx.o fragLength

%.o: %.cpp #../lib/libgab/utils.o  
	${CXX} ${CXXFLAGS} $^ -o $@

fragLength:	fragLength.o   ../lib/libgab/utils.o   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

../lib/htslib/libhts.a:  ../lib/htslib/hts_internal.h
	cd ../lib/htslib/ && make && cd ../..

../lib/htslib/hts_internal.h:
	rm -rfv ../lib/htslib/
	git clone --recursive https://github.com/samtools/htslib.git ../lib/htslib/

../lib/samtools/bedidx.o:  ../lib/samtools/sam.h
	cd ../lib/samtools/ && make && cd ../..

../lib/samtools/sam.h:
	rm -rf ../lib/samtools/
	git clone --recursive https://github.com/samtools/samtools.git ../lib/samtools/

../lib/libgab/utils.o:  ../lib/libgab/utils.h
	cd ../lib/libgab/ &&  make utils.o && make -C gzstream/ && cd ../..

../lib/libgab/utils.h:
	rm -rfv ../lib/libgab/
	mkdir -p ../lib/libgab/
	git clone --depth 1 https://github.com/grenaud/libgab.git ../lib/libgab/


clean :
	rm -f *.o bam2prof diffprof


