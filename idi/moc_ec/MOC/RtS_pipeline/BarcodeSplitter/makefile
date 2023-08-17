# compiler to use
CC = c++
CFLAGS=-O3 -std=c++14  
OFLAGS=-lpthread -fopenmp

LOCALPATH=/broad/IDP-Dx_work/nirmalya/local
INC = -I${LOCALPATH}/include 

BOOSTLIBS=${LOCALPATH}/lib/libboost_serialization.a ${LOCALPATH}/lib/libboost_regex.a -lboost_iostreams -lz
PROG_OPT_LIB=${LOCALPATH}/lib/libboost_program_options.a
STXXLLIBS=${LOCALPATH}/lib/libstxxl.a
 
all: clean tools
	
tools:
	#$(CC) $(CFLAGS) $(INC) dict_builder.cpp -o dict_builder $(BOOSTLIBS) $(PROG_OPT_LIB)
	#$(CC) $(CFLAGS) $(INC) index_splitter.cpp -o index_splitter $(BOOSTLIBS) $(PROG_OPT_LIB)
	#$(CC) $(CFLAGS) $(INC) dict_builder_test.cpp -o dict_builder_test $(BOOSTLIBS) $(PROG_OPT_LIB)
	#$(CC) $(CFLAGS) $(INC) barcode_splitter.cpp -o bc_splitter $(BOOSTLIBS) $(PROG_OPT_LIB)
	$(CC) $(CFLAGS) $(INC) barcode_splitter_rts.cpp -o bc_splitter_rts $(BOOSTLIBS) $(PROG_OPT_LIB)
	#$(CC) $(CFLAGS) $(INC) barcode_splitter_rts_se.cpp -o bc_splitter_rts_se $(BOOSTLIBS) $(PROG_OPT_LIB)
	#$(CC) $(CFLAGS) $(INC) fastq_gz_demo.cpp -o fastq_gz_demo $(BOOSTLIBS) $(PROG_OPT_LIB)
	#$(CC) $(CFLAGS) $(INC) boostgz.cc -o boostgz $(BOOSTLIBS) $(PROG_OPT_LIB)
	
clean:
	rm -f bkLoad bkSearch
	rm -Rf bkSearch.dSym bkLoad.dSYM
