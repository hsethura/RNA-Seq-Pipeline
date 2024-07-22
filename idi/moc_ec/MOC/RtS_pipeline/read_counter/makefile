# compiler to use
CC = c++
CFLAGS=-O3 -std=c++14 -g

LOCALPATH=/broad/IDP-Dx_work/nirmalya/local
INC = -I${LOCALPATH}/include 
LIBDIR=${LOCALPATH}/lib/
#LIBS=${LIBDIR}/libhts.so  -lz -lboost_program_options -lboost_regex
LIBS=-lhts  -lz -lboost_program_options -lboost_regex

 
all: tools
	
tools:
	#$(CC) $(CFLAGS) $(INC) sam_fragcount.cpp -o sam_fragcount  $(LIBS) 
	#$(CC) $(CFLAGS) $(INC) frag_to_gene_count.cpp -o frag_to_gene_count  $(LIBS)
	#$(CC) $(CFLAGS) $(INC) metrics_gen.cpp -o metrics_gen  $(LIBS) 
	#$(CC) $(CFLAGS) $(INC) filter_sam.cpp -o filter_sam  $(LIBS) 
	#$(CC) $(CFLAGS) $(INC) paired_only.cpp -o paired_only  $(LIBS) 
	#$(CC) $(CFLAGS) $(INC) remove_dup.cpp -o remove_dup  $(LIBS) 
	$(CC) $(CFLAGS) $(INC) remove_bu.cpp -o remove_bu  $(LIBS) 
	
