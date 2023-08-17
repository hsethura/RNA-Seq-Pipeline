#!/bin/sh

PROJ_ID=$1
shift

USER_ID=`id | cut -d"(" -f2 | cut -d")" -f1`


### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

SCRIPT_OPTIONS=$@

CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $SCRIPT_OPTIONS`
TOP_GENE_NUM=`extract_option -top_gene 500 1 $SCRIPT_OPTIONS`
MIN_LEN=`extract_option -min_len 800 1 $SCRIPT_OPTIONS`
ADD_5=`extract_option -add_5 20 1 $SCRIPT_OPTIONS`
ADD_3=`extract_option -add_3 30 1 $SCRIPT_OPTIONS`

### set paths to directories, files, and scripts from config file

RESULTS_PATH=`config_read $CONFIG_FILE Results_path`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`

INDEX1_BARCODES=`config_read $CONFIG_FILE P7_barcodes`
SCRIPT_DIR=`config_read $CONFIG_FILE pipe_script_dir`

COORD_DIR=$TEMP_PATH$USER_ID"/"$PROJ_ID"/patho_result/"$PROJ_ID"/temp_bamdir/"
GFF_DIR=$TEMP_PATH$USER_ID"/"$PROJ_ID"/datadir/"

TEMP_DIR=$TEMP_PATH$USER_ID"/"

mkdir -p $TEMP_DIR

ALL_COORDS=`ls -lrt $COORD_DIR* | grep coords.txt | awk '{print $9}'`

ALL_ACCS=`for m in $ALL_COORDS
		do
			cat $m | awk '{print $1}' | sort | uniq
			exit
		done`


for COORDS in $ALL_COORDS
do
	for ACC in $ALL_ACCS
	do
		
		ROOT=`echo $COORDS | rev | cut -d"/" -f1 | rev | sed 's/_coords.txt//g' `
		CP_OUT=$TEMP_DIR$ROOT"_CP.txt"
		DIST_TEMP=$TEMP_DIR$ROOT"_DIST.txt"
		DIST_TEMP=$TEMP_DIR$ROOT"_DIST.txt"
		
		GFF_FILE=$GFF_DIR$ACC"_GENES.gff"
		GFF_PARSED=$TEMP_DIR$ACC"_GENES_parsed.gff"
		
		cat $GFF_FILE | awk '{
									y=split($9, ar, ";")
									for(i=1; i < y+1; i++)
									{	
										if(ar[i] ~ /full_tag=/)
										{	
											split(ar[i], pr, "=")
											name=pr[2]
										}
									}
								print name, $4, $5, $7
							
							}' > $GFF_PARSED
									
		$SCRIPT_DIR"/shell_scripts/COORDS_TO_CP" $COORDS $ACC Y > $CP_OUT
 		$SCRIPT_DIR"/shell_scripts/cov_dist" $GFF_PARSED $CP_OUT 20 $ADD_5 $ADD_3 | awk -v MIN_LEN=$MIN_LEN '{y=split($1, ar, ":"); len=ar[y]; if(len>=MIN_LEN) print $0}' | sort -k2n | tail -$TOP_GENE_NUM > $DIST_TEMP
		
		printf "%s\t", $ROOT
		
		cat $DIST_TEMP | awk '{
				for(i=3; i < NF+1; i++)
				{
					if(NR==1)
						total[i]=0
					total[i]=$2*($i/100)
					printf "%s\t", total[i]
				}
			print ""
		
		}' | tail -1  

		echo $DIST_TEMP

	done
done
