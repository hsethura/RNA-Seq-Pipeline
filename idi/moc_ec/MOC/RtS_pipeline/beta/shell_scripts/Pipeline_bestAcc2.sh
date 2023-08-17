#!/bin/sh

KEY_FILE=$1
shift
BAM_DIR=$1"/"
shift
TEMP_PATH=$1"/"
shift
CONFIG_FILE=$1
shift
SCRIPT_DIR=$1"/"
shift

source "/broad/software/scripts/useuse"
reuse -q UGER
reuse -q Java-1.6
reuse -q BWA
reuse -q BLAST+

##############  FUNCTIONS  #################

########## function for extracting field numbers from key headers
 
find_field()
{
	INPUT_FILE=$1
	field_title=$2
	req_field=$3
	
	cat $INPUT_FILE | grep -v "#" | head -1 | awk -F"\t" -v f=$field_title -v req=$req_field '{
								x=0
								for(i=1; i < NF + 1; i++)
									if($i == f)
										x=i
								if(x==0 && req=="Y")
									print f";exit"
								if(x!=0 || (x==0 && req=="N"))
									print x
								
							}'
}

########## function for pulling out paths, header names from config file 

config_read ()
{

	CONFIG_FILE=$1
	VAR_NAME=$2
	
	F=`echo "0" && cat $CONFIG_FILE | grep -w $VAR_NAME":" | awk '{ print $2 }'` 

	N=`echo $F | awk '{print $NF}'`
	
	if [ $N == 0 ];then
		echo $VAR_NAME" not found in "$CONFIG_FILE
		echo $F
		exit
	fi
	echo $N
}
########## function for pulling out paths, header names from config file 

config_find_field ()
{

	CONFIG_FILE=$1
	INPUT_FILE=$2
	VAR_NAME=$3
	req_field=$4

	F=`echo "0" && cat $CONFIG_FILE | grep -w $VAR_NAME":" | awk '{ print $2 }'` 

	N=`echo $F | awk '{print $NF}'`
			
	cat $INPUT_FILE | grep -v "#" | head -1 | awk -F"\t" -v f=$N -v req=$req_field '{
							x=0
							for(i=1; i < NF + 1; i++)
								if($i == f)
									x=i
							if(x==0 && req=="Y")
								print f";exit"
							if(x!=0 || (x==0 && req=="N"))
								print x
							
						}'

}

############################  SGE_test ###############################

SGE_test ()
{	
	SGE_OUT=$1
	
	jobs_submit=`cat $SGE_OUT`
	
	typeset -i NUM
	
	while :
	do
			
		job_running=`qstat | awk '{if(NR > 2 ) print $1}'`
		NUM=`echo $job_running $jobs_submit | tr ' ' '\n' | sort | uniq -c | awk '{if($1==2)print $0}' | wc -l`
		if [ $NUM == 0 ];then
		
			echo "qstat jobs are complete"
			break
		fi
		
	done
}

########################################################################

########  function for extracting options from command line
 
		extract_option ()
		{
			option_name=$1 
			shift
			default=$1
			shift
			
			echo $@ | awk -v def=$default -v name=$option_name '{
		
				for(i=1; i < NF+1; i++)
				{
					if(i==1)			# intitalize to $default in case no option found
						print def
					if($i==name)
						print $(i+1)
				}
		
			}' | tail -1 
			
			
		}

####################################################

USER=`id | awk '{print $1}' | cut -d"(" -f2 | sed 's/)//g'` 

options=$@" -"

STRAIN=`extract_option -STRAIN - $options`
USERID=`extract_option -USERID $USER $options`

############  Extract field IDs from key headers  #################
 
ls -lrt $CONFIG_FILE

NCBI_PATH=`config_read $CONFIG_FILE patho_dbpath`
ALN_DIR=`config_read $CONFIG_FILE Align_path`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
#SCRIPT_DIR=`config_read $CONFIG_FILE Script_dir`
KEY_DIR=`config_read $CONFIG_FILE Key_base`


TEMP_DIR=$TEMP_PATH"/"$USERID"/best_acc/"

mkdir -p $TEMP_DIR

echo $TEMP_DIR 

ls -lrt $KEY_FILE

INDEX_F=`config_find_field $CONFIG_FILE $KEY_FILE P7 N`
GENUS_F=`config_find_field $CONFIG_FILE $KEY_FILE Genus Y`
STRAIN_F=`config_find_field $CONFIG_FILE $KEY_FILE Strain N`
BC_F=`config_find_field $CONFIG_FILE $KEY_FILE bc N` 
PROJID_F=`config_find_field $CONFIG_FILE $KEY_FILE Proj N`
SAMPID_F=`config_find_field $CONFIG_FILE $KEY_FILE ID N`
REF_F=`config_find_field $CONFIG_FILE $KEY_FILE Ref_accession N`
SEQFILE_F=`config_find_field $CONFIG_FILE $KEY_FILE path_to_seq N`

echo $INDEX_F
echo $GENUS_F
echo $STRAIN_F
echo config_find_field $CONFIG_FILE $KEY_FILE Genus Y
exit

p=`echo $INDEX_F $GENUS_F $STRAIN_F $BC_F $PROJID_F $SAMPID_F $PROJID_F $REF_F`

for m in $p
do
	
	echo $m
	if [ $m == "0;exit" ];then
		echo "Problem with one of the fields"
		echo "INDEX_F" $INDEX_F 
		echo "GENUS_F" $GENUS_F 
		echo "STRAIN_F" $STRAIN_F 
		echo "BC_F" $BC_F 
		echo "PROJID_F" $PROJID_F 
		echo "SAMPID_F" $SAMPID_F
		echo "PROJID_F" $PROJID_F 
		echo "REF_F" $REF_F
		echo "SEQFILE_F" $SEQFILE_F
		exit 
	fi
done

if [ $STRAIN == "-" ];then
	ALL_STRAIN_FS=`cat $KEY_FILE |  awk -F"\t" '{if($0 !~ /\#/ && $0 !~ /Project_ID/ && $0 !~ /mixed/) print $'$STRAIN_F'}' | sort | uniq` 
else
	ALL_STRAINS=`echo $STRAIN | sed 's/:/ /g'`
	ALL_STRAIN_FS=`for m in $ALL_STRAINS
		do
			cat $KEY_FILE |  awk -F"\t" -v STRAIN_F=$STRAIN_F '{if($0 !~ /\#/ && $0 !~ /Project_ID/ && $0 !~ /mixed/) print $STRAIN_F}' | grep $m  
		done | sort | uniq`
fi

echo "All Strains: "$ALL_STRAIN_FS


SPLIT_FILES=""


for m in $ALL_STRAIN_FS
do
	FILES=`cat $KEY_FILE | grep $m | head -1 | grep -v "#" | awk -F"\t" '{
	
				split($'$REF_F', ar, ";")
				REF=ar[1]
				
				FILES="'$BAM_DIR'/"$'$SAMPID_F'"_R2.fastq"
				#FILES="'$BAM_DIR'/"$'$PROJID_F'"/"$'$SAMPID_F'"_R2.fastq"
				print REF";"$'$GENUS_F'";"$'$STRAIN_F'";"FILES				
				#print $'$BC_F'
		
			}'` 
	
	SPLIT_FILES=$SPLIT_FILES" "$FILES
done	 

echo $SPLIT_FILES


for m in $SPLIT_FILES
do
	echo $m

	REF_F=`echo $m | cut -d";" -f1`
	if [ $REF_F == "Unknown" ] || [ $REF_F == "Unknown;" ] || [ $STRAIN != "-" ];then
		GENUS_F=`echo $m | cut -d";" -f2`
		STRAIN_F=`echo $m | cut -d";" -f3` 
		SEQ_FILE=`echo $m | cut -d";" -f4`		
		FILE_TYPE=`echo $SEQ_FILE | rev | cut -d"." -f1 | rev`
		
		echo $SEQ_FILE 
			
		if [ $FILE_TYPE == "bam" ];then
			
			ROOT=`echo $SEQ_FILE | awk '{y=split($1, ar, "/"); print ar[y]}' | sed -e 's/.gz//g' -e 's/.fastq//g' -e 's/.bam//g'`
			if [ ! -s $TEMP_DIR"/"$ROOT".fastq.1" ];then
				sh $SCRIPT_DIR"ALN_BAM_to_FASTQ.sh" $SEQ_FILE $ROOT 1 1 $TEMP_DIR
			fi
			SEQ_FILE=$TEMP_DIR"/"$ROOT".fastq.1"
		fi

		echo "sh $SCRIPT_DIR"Pipeline_align_test2.sh" $STRAIN_F $GENUS_F 10000 $SEQ_FILE $ALN_DIR $TEMP_DIR $NCBI_PATH $SCRIPT_DIR"
		sh $SCRIPT_DIR"Pipeline_align_test2.sh" $STRAIN_F $GENUS_F 10000 $SEQ_FILE $ALN_DIR $TEMP_DIR $NCBI_PATH $SCRIPT_DIR
	
	fi
done

 
