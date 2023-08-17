#!/bin/sh

DIR=$1


MOC_ID="none"

source /broad/software/scripts/useuse

use Python-2.7

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

rev_comp ()
{
	
	FILE=$1
	
	cat $FILE | awk -v STRAND=$STRAND '{ 
				
					
							y=length($2)
							printf "%s\t", $1
						
							for(t=y; t > 0; t--)
							{
								p=substr($2, t, 1)
						
								if(p=="A")
									printf "T"
								if(p=="C")
									printf "G"
								if(p=="G")
									printf "C"
								if(p=="T")
									printf "A"
								if(p=="a")
									printf "t"
								if(p=="c")
									printf "g"
								if(p=="g")
									printf "c"
								if(p=="t")
									printf "a"
							}
							print ""
		}'
}

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`
OW_SYM=`extract_option -ow_sym N 1 $@`

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

mkdir -p MOCSDB_DIR

### set paths to directories, files, and scripts from config file

TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
GLOCAL_PATH=`config_read $CONFIG_FILE gdrivelocal_path`
GDRIVE_SCRIPT=`config_read $CONFIG_FILE gdrive_script`
TEMP_DIR=~/HPV/
mkdir -p $TEMP_DIR


DATA_DIR=`extract_option -data_dir $PARSE_DIR 1 $@`
PARSE=`extract_option -parse Y 1 $@`
IMPORT_GS=`extract_option -import_gs Y 1 $@`

GMOCPM_PATH=`config_read $CONFIG_FILE gdrivepm_path`
GWB_DIR=$GMOCPM_PATH"/Templates_and_workbooks/RtS/IDing_and_submission_WBs/"
LOCAL_DIR=$GLOCAL_PATH"/"$GWB_DIR

RAW_SYM_PATH=`extract_option -symlink $RAWSYM_PATH 1 $@`
RAW_SEQ_PATH=`extract_option -raw_seq_path $RAWSYM_PATH 1 $@`

echo "RAW_SEQ_PATH: "$RAW_SEQ_PATH
echo "$RAW_SYM_PATH: "$RAW_SYM_PATH

mkdir -p $LOCAL_DIR

echo "cd $GLOCAL_PATH"
cd $GLOCAL_PATH
pwd 


mkdir -p $BCSDB_DIR
echo "$SCRIPTS_DIR"GS_import.py" -s $BCS_GID -t "NEW barcodes and indexes" -p $BCS_FILE_NAME --Key_dir $BCSDB_DIR"
$SCRIPTS_DIR"GS_import.py" -s $BCS_GID -t "NEW barcodes and indexes" -p $BCS_FILE_NAME --Key_dir $BCSDB_DIR -S $BCS_FILE_SUFF

BCS_FILE=$BCSDB_DIR"/"$BCS_FILE_NAME$BCS_FILE_SUFF
ls -lrt $BCS_FILE


echo "FIELD_HEADER $BCS_FILE HPV_Combined"
BC_F=`FIELD_HEADER $BCS_FILE HPV_Combined`

ALL_BCSEQ=`cat $BCS_FILE | sed 1d | awk -F"\t" '{print $'$BC_F'}' | sort | uniq`
echo $ALL_BCSEQ

ALL_FASTQS=`ls -lrt $DIR | grep fastq.gz | awk '{print $NF}'`
echo $ALL_FASTQS

ALL_IDS=`ls -lrt $DIR | grep fastq.gz | awk '{print $NF}' | sed 's/_R[0-9]_00[0-9].fastq.gz//g' | sort | uniq `
echo $ALL_IDS

ALL_MET_FILES=""

for ID in $ALL_IDS
do
	echo $ID
	FQ1=`ls -lrt $DIR"/"$ID*"R1"* | awk '{print $NF}'`
	FQ2=`ls -lrt $DIR"/"$ID*"R2"* | awk '{print $NF}'`
	FQ1_TEMP=$TEMP_DIR"/"$ID"_R1_temp.fastq"
	FQ2_TEMP=$TEMP_DIR"/"$ID"_R2_temp.fastq"
	FQ_TEMP=$TEMP_DIR"/"$ID"_COMBO_temp.fastq"
	

	zcat $FQ1 | awk '{if(NR % 4==2)print substr($1,1,100)}' > $FQ1_TEMP
	zcat $FQ2 | awk '{if(NR % 4==2)print substr($1,1,100)}' > $FQ2_TEMP
	paste $FQ1_TEMP $FQ2_TEMP > $FQ_TEMP
	
	OUT_METS=$TEMP_DIR"/"$ID"_metrics.txt"
	ALL_MET_FILES=$ALL_MET_FILES" "$OUT_METS
		
	for BC in $ALL_BCSEQ
	do
		WID=`echo $BC | sed 's/_/ /g' | awk '{print $1}'`
		BC1=`echo $BC | sed 's/_/ /g' | awk '{print $2}'`
		BC2=`echo $BC | sed 's/_/ /g' | awk '{print $3}'`
	
		BC1_LEN=`echo $BC1 | awk '{print length($1)}'`
		BC2_LEN=`echo $BC2 | awk '{print length($2)}'`

		OUT_READS=$TEMP_DIR"/"$ID"_"$WID".txt"

		echo "" | sed 1d > $OUT_READS
		
		cat $FQ_TEMP | awk -v total=0 -v BC1_LEN=$BC1_LEN -v BC2_LEN=$BC2_LEN -v BC1=$BC1 -v BC2=$BC2 '{if($1 ~ BC1 && $2 ~ BC2) print substr($1,BC1_LEN+25,1000), substr($2,BC2_LEN+25,1000)}' | rev_comp  >> $OUT_READS

		ls -lrt $OUT_READS
		
		READ_NUM=`cat $OUT_READS | wc -l`
		NUM_HPV16_1=`cat $OUT_READS | awk '{print $1}' | grep -e "TGTCATTATGTGCTGCCATATCTACTTCAGAA" -e "TATTTGTTACTGTTGTTGATACTACACGCA"| wc -l`
		PCNT_HPV16_1=`echo $NUM_HPV16_1 $READ_NUM | awk '{print $1/$2*100}'`
		
		NUM_HPV16_2=`cat $OUT_READS | awk '{print $1}' | grep -e "GTACAAATATGTCATTATGTGCTGCCATATCTACTTCAGAA" -e "TATTTGTTACTGTTGTTGATACTACACGCA" | wc -l`
		PCNT_HPV16_2=`echo $NUM_HPV16_2 $READ_NUM | awk '{print $1/$2*100}'`

# 		PCNT_TOP2=`cat $OUT_READS | awk '{print $2}' | sort -T /idi/moc_ec | uniq -c | sort | tail -1 | awk '{print $1/'$READ_NUM'*100}'`
# 		echo $ID $WID $READ_NUM $PCNT_TOP1 $PCNT_TOP2 >> $OUT_METS
		echo $ID $WID $READ_NUM $PCNT_HPV16_1 >> $OUT_METS

	done
	ls -lrt $OUT_METS

done

BATCH_ID=`basename $DIR`
FINAL_READ_MET=$TEMP_DIR$BATCH_ID"_read_metrics.txt"

echo "WellID" > $FINAL_READ_MET
paste $ALL_MET_FILES | awk '{print $2}' >> $FINAL_READ_MET

for MET_FILE in $ALL_MET_FILES
do
	echo $MET_FILE
	ID_ID=`cat $MET_FILE | awk '{print $1}' | head -1`
	echo "/t"$ID_ID > temp.txt
	cat $MET_FILE | awk '{print $3}' >> temp.txt 
	cat $FINAL_READ_MET > temp2.txt
	paste temp2.txt  temp.txt >> $FINAL_READ_MET
done

ls -lrt $FINAL_READ_MET
echo $FINAL_READ_MET




ls -lrt ~/combo_fastq.txt
### change permissions for Results and temp dirs
change_perms $LOCAL_DIR $MOVE_DB $BCSDB_DIR $SYM_DIR $MOCSDB_DIR
