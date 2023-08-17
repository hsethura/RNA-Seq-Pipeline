# compiler to use
CC = c++
CFLAGS=-O3 -std=c++14 -g

LOCALPATH=/broad/IDP-Dx_work/nirmalya/local
INC = -I${LOCALPATH}/include
LIBDIR=${LOCALPATH}/lib/
#LIBS=${LIBDIR}/libhts.so  -lz -lboost_program_options -lboost_regex
LIBS=-lhts  -lz -lboost_program_options -lboost_regex -lboost_iostreams -lz


all: tools
    
tools:
	$(CC) $(CFLAGS) $(INC) sample_fastq.cpp -o sample_fastq  $(LIBS) 
	$(CC) $(CFLAGS) $(INC) seed_fgen.cpp -o seed_fgen  $(LIBS) 

