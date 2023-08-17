# compiler to use
CC = c++
CFLAGS=-O3 -std=c++14

LOCALPATH=/broad/IDP-Dx_work/nirmalya/local
INC = -I${LOCALPATH}/include
LIBDIR=${LOCALPATH}/lib/
LIBS=-lz -lboost_program_options -lboost_iostreams

 
all: clean tools
	
tools:
	$(CC) $(CFLAGS) $(INC) read_trimmer.cpp -o read_trimmer $(LIBS)
	
clean:
