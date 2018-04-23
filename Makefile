
# include SeqAn libraries
#CXXFLAGS+=-I/nfs/prog/bioinfo/u/arnaos/seqanProject/seqan-1.4.1/core/include
MAX_KMER_SIZE=64
CXXFLAGS+=-DMAX_KMER_SIZE=$(MAX_KMER_SIZE) -I../libraries/seqan-1.4.2/include

# RELEASE build 
CXXFLAGS+= -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1 -g
LDLIBS=-lz -std=c++0x

# set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x

all: bamqualcheck

OBJECTS = lsb.o RepHash.o Kmer.o KmerIterator.o hash.o

bamqualcheck: bamqualcheck.o $(OBJECTS)
	$(CXX) $(OBJECTS) bamqualcheck.o $(LDFLAGS) -o bamqualcheck $(LDLIBS)


lsb.o: lsb.cpp lsb.hpp
bamqualcheck.o: bamqualcheck.cpp StreamCounter.hpp
Kmer.o: Kmer.cpp Kmer.hpp
KmerIterator.o: KmerIterator.cpp KmerIterator.hpp
hash.o: hash.cpp hash.hpp
RepHash.o: RepHash.cpp RepHash.hpp

clean:
	rm -f *.o ./bamqualcheck
