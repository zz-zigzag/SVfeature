
CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

SEQAN = ${NGS_TOOL}/seqan-library-2.0.0
INC = -I${SEQAN}/include

CFLAGS = -fexceptions -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -DSEQAN_HAS_ZLIB=1
RESINC = 
LIBDIR = 
LIB = -lz -lpthread
LDFLAGS = -static

INC_RELEASE = $(INC)
CFLAGS_RELEASE = $(CFLAGS) -O2 -DSEQAN_ENABLE_DEBUG=0
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -s

all: release
release: main.o
	$(LD) $(LIBDIR_RELEASE) -o SVfeature main.o $(LDFLAGS_RELEASE) $(LIB_RELEASE)

main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o main.o

clean: 
	rm main.o SVfeature

.PHONY: clean

