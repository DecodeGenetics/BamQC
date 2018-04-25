TARGET = bamqualcheck
BUILD_DIR = ./build
SRC_DIR = ./src

HDRS := $(wildcard $(SRC_DIR)/*.hpp) $(wildcard $(SRC_DIR)/*.h)
SRCS := $(wildcard $(SRC_DIR)/kmerstream/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SRCS:.cpp=.o))

# Set k-mer size for k-mer stream
MAX_KMER_SIZE=64
CXXFLAGS+=-DMAX_KMER_SIZE=$(MAX_KMER_SIZE)

# Include SeqAn libraries
CXXFLAGS+=-I../libraries/seqan-1.4.2/include

# Set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x

# RELEASE build 
CXXFLAGS+= -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1 -g
LDLIBS=-lz -std=c++0x

all: $(TARGET)

$(TARGET): $(BUILD_DIR) $(BUILD_DIR)/$(TARGET).o $(OBJS)
	$(CXX) $(BUILD_DIR)/$(TARGET).o $(OBJS) -o $@ $(LDLIBS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR) $(BUILD_DIR)/kmerstream

$(BUILD_DIR)/$(TARGET).o: $(SRC_DIR)/$(TARGET).cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/kmerstream/RepHash.o: $(SRC_DIR)/kmerstream/RepHash.cpp $(SRC_DIR)/kmerstream/RepHash.hpp $(SRC_DIR)/kmerstream/mersennetwister.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(BUILD_DIR)/$(TARGET).o $(OBJS) $(TARGET)
	rmdir build/kmerstream build
