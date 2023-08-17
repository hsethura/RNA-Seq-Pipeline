#!/bin/sh

MOC_ID=$1
PROJ_ID=$2
READ_DIR=$3
OUT_DIR=$4
ACC=$5
COND1_ID=$6
COND2_ID=$7
COND1_TAG=$8
COND2_TAG=$9


read_file_maker ()
{
	READ_FILE=$1
	shift
	TEMP_DIR=$1
	shift
	READ_OUTPUT=$1
	shift

	ALL_COND=$@
	
	TEMP_OUTFILE1=$TEMP_DIR"/temp1.txt"
	TEMP_OUTFILE2=$TEMP_DIR"/temp2.txt"
	TEMP_OUTFILE3=$TEMP_DIR"/temp3.txt"

	cat $READ_FILE | awk -F"\t" '{print $1}' > $TEMP_OUTFILE1

	for m in $ALL_COND
	do
		NUM_MATCH=`cat $READ_FILE | head -1 | grep -w $m  | wc -l`
		echo $NUM_MATCH, $m, $READ_FILE 
		if [ $NUM_MATCH -eq 0 ];then
	
			echo "No match found for "$m" in "$READ_FILE"!"
			exit
		fi
	
	
	
		cat $READ_FILE | awk -F"\t" -v name=$m -v TEMP_OUTFILE2=$TEMP_OUTFILE2 '{
													
														if(NR == 1)
															for(i=1; i < NF+1; i++)
																if($i == name)
																	y=i

														print $y > TEMP_OUTFILE2
 
																											
													}'
	
	
		cat $TEMP_OUTFILE1 > $TEMP_OUTFILE3
	
		paste $TEMP_OUTFILE3 $TEMP_OUTFILE2 > $TEMP_OUTFILE1
	
		echo $TEMP_OUTFILE1, $TEMP_OUTFILE2
	done 


	cat $TEMP_OUTFILE1 > $READ_OUTPUT
}

###############  function for extracting options ###############

# Extracting options from the command line
	extract_option () 
	{
		option_name=$1
		shift
		default=$1
		shift 
		max=$1
		shift 
	
		typeset -i num
		num=`echo $@ | wc -w`
		if [ $num -eq 0 ];then
		
			echo $default
			exit
			
		fi
		
		echo $@ | awk -v def=$default -v name=$option_name -v flag="N" '{
	
			for(i=1; i < NF+1; i++)
			{
				if($i==name)
				{
					flag="Y"
					for(n=1; n < (NF+1)-i; n++)
					{
						y=substr($(n+i),1,1)
						##// if next option encountered, exit
						if(y=="-")
							exit
						print $(i+n)
					}
				}
			
				if(flag=="N" && i==NF)
					print def
			
			}
	
		}' | tail -$max
	}
###########################################################

# Finding field #'s using field headers

find_field()
{
	INPUT_FILE=$1
	field_title=$2
	req_field=$3
	
	cat $INPUT_FILE | head -1 | awk -F"\t" -v f=$field_title -v req=$req_field '{
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
#################################################

# Making DESeq2 R files

DESeq2R_maker ()
{
	READ_OUTPUT=$1
	KEY_OUTPUT=$2
	DS2_RESULTS=$3
	DS2_MA=$4
	tag=$5
	test=$6

	echo "library(DESeq2)"
	echo "counts = read.table(\""$READ_OUTPUT"\", header=T, row.names = 1, sep = \"\t\", quote = \"\", comment.char = \"\")" 
	echo "key = read.table(\""$KEY_OUTPUT"\", header=T)" 
	echo "rownames(key) = key[,1]"
	echo "des = DESeqDataSetFromMatrix(counts, key, ~Group)" 
	
	if [ $test == "WALD" ];then
		echo "des = DESeq(des, modelMatrixType = \"expanded\")" 
	fi
	
	if [ $test == "LRT" ];then
		echo "des = DESeq(des, test=\"LRT\", reduced=~1)"
	fi

	echo "res = results(des, contrast = c(\"Group\", \"A\", \"B\"))" 
	echo "write.table(res, file = \""$DS2_RESULTS"\",  sep = \"\t\", quote = F, row.names = T)" 
	echo "x = order(res[,2], decreasing = TRUE)"
	echo "max = res[x[1],2]"
	echo "x = order(res[,2], decreasing = FALSE)"
	echo "min = res[x[1],2]"
	echo "pdf(\""$DS2_MA"\", width=4, height=4)"
	echo "plotMA(res, main=\""$tag"\", ylim=c(min,max))"
	echo "dev.off()"

}

####################################

# Making DESeq1 R files

DESeq1R_maker ()
{
	READ_OUTPUT=$1
	KEY_OUTPUT=$2
	DS1_RESULTS=$3
	DS2_MA_CDS=$4
	tag=$5
	numA=$6
	numB=$7
	
	echo "library(DESeq)" 
	echo "counts <- read.delim(\""$READ_OUTPUT"\", row.names=1)"
	printf "conditions <- c("
	echo "" | awk -v n=$numB '{
						for(i=1; i < n+1; i++)
							printf "\"B\","
					}' 
	echo "" | awk -v n=$numA '{
						for(i=1; i < n; i++)
							printf "\"A\","
						printf "\"A\""
					}' 

	echo ")"

	echo "countsnew1 <- newCountDataSet(counts, conditions)"
	echo "countsnew2 <- estimateSizeFactors(countsnew1)"
	echo "countsnew3 <- estimateDispersions(countsnew2)"
	echo "nbintest <- nbinomTest(countsnew3, \"A\", \"B\")"
	echo "write.table(nbintest, file=paste(\""$DS1_RESULTS"\", sep=\"""\"),"
	echo "quote=FALSE, sep=\"\t\", row.names=FALSE)"

}


##################################################

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

source "/broad/software/scripts/useuse"
reuse -q UGER
reuse -q Java-1.8
reuse -q R-3.2

CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $@`

# Get paths to dirs scripts from config file
read_config $CONFIG_FILE 

### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.
path_suff $@

# MAKING READ FILE

min_fold=`extract_option -f 2 1 $@`
max_P=`extract_option -P 0.05 1 $@`
include=`extract_option -i - 100 $@`
compare=`extract_option -c - 2 $@`
omit=`extract_option -o ABCDEFGHIJKLMNOP 100 $@`
tag=`extract_option -tag - 1 $@`
test=`extract_option -test WALD 1 $@`


CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $SCRIPT_OPTIONS`

# Get paths to dirs scripts from config file
read_config $CONFIG_FILE 

KEY_FILE=$KEY_DIR"/"$MOC_ID"_key.txt"
TEMP_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/DE/"

echo $KEY_FILE
echo $MOC_ID

TEMP_KEY_FILE=$TEMP_DIR"/"$MOC_ID"_temp_key.txt"

TAG=$COND1_TAG"_vs_"$COND2_TAG"_"$test
DS_DIR=$TEMP_DIR"/temp_files/"

mkdir -p $DS_DIR 
mkdir -p $OUT_DIR
mkdir -p $TEMP_DIR

READ_FILE_TEMP=$READ_DIR"/"$PROJ_ID"_"$ACC"_counts.tsv"
READ_FILE_CDS=$TEMP_DIR$PROJ_ID"=REP_Genes_"$ACC"_CDS_Reads.txt"
READ_FILE_ALL=$TEMP_DIR$PROJ_ID"=REP_Genes_"$ACC"_AllFeatures_Reads.txt"
DS2_RFILE_CDS=$TEMP_DIR$PROJ_ID"-"$TAG"_DS2_CDS.R"
DS2_RFILE_ALL=$TEMP_DIR$PROJ_ID"-"$TAG"_DS2_AllFeatures.R"

if [ -s $READ_FILE_TEMP ];then
	cat $READ_FILE_TEMP | grep -v tRNA: | grep -v IGR: | grep -v AS_ | grep -v rRNA: > $READ_FILE_CDS
	cat $READ_FILE_TEMP > $READ_FILE_ALL
else
	echo $READ_FILE_TEMP" not found"
fi
 
ls -lrt $KEY_FILE
ls -lrt $READ_FILE_CDS


if [ ! -s $READ_FILE_CDS ];then

	echo $READ_FILE_CDS" not found!"

fi
if [ ! -s $KEY_FILE ];then

	echo $KEY_FILE" not found!"

fi
			
ALL_COND=`echo $COND1_ID $COND2_ID | sed 's/,/ /g'`

echo $ALL_COND

KEY_OUTPUT=$DS_DIR$PROJ_ID"-"$TAG"_key.txt"
READ_OUTPUT_CDS=$DS_DIR$PROJ_ID"-"$TAG"_reads_CDS.txt"
READ_OUTPUT_ALL=$DS_DIR$PROJ_ID"-"$TAG"_reads_AllFeatures.txt"

read_file_maker $READ_FILE_CDS $TEMP_DIR $READ_OUTPUT_CDS $ALL_COND
read_file_maker $READ_FILE_ALL $TEMP_DIR $READ_OUTPUT_ALL $ALL_COND

echo $READ_OUTPUT_CDS
echo  $READ_OUTPUT_ALL

echo "Sample_ID	Condition	Group" > $KEY_OUTPUT

f=`echo $COND1_ID | sed 's/,/ /g'`

for m in $f
do
	echo $m"	"$COND1_TAG"	A" >> $KEY_OUTPUT
done

f=`echo $COND2_ID | sed 's/,/ /g'`
for m in $f
do
	echo $m"	"$COND2_TAG"	B" >> $KEY_OUTPUT
done

echo $READ_OUTPUT
echo $KEY_OUTPUT

COND1_NUM=`echo $COND1_ID | sed 's/,/ /g' | wc -w`
COND2_NUM=`echo $COND2_ID | sed 's/,/ /g' | wc -w`
OUTPUT_TAG=$TAG"_"$COND1_NUM","$COND2_NUM

tag2="_"$numA"vs"$numB
DS2_RESULTS_CDS=$TEMP_DIR"DESeq2_results_"$PROJ_ID"-"$OUTPUT_TAG"_CDS.txt"
DS2_OUTPUT_CDS=$TEMP_DIR"DESeq2_out_"$PROJ_ID"-"$OUTPUT_TAG"_CDS.txt"
DS2_MA_CDS=$TEMP_DIR"DESeq2_"$PROJ_ID"-"$OUTPUT_TAG"_MA_CDS.pdf"

DS2_RESULTS_ALL=$TEMP_DIR"DESeq2_results_"$PROJ_ID"-"$OUTPUT_TAG"_AllFeatures.txt"
DS2_OUTPUT_ALL=$TEMP_DIR"DESeq2_out_"$PROJ_ID"-"$OUTPUT_TAG"_AllFeatures.txt"
DS2_MA_ALL=$TEMP_DIR"DESeq2_"$PROJ_ID"-"$OUTPUT_TAG"_MA_AllFeatures.pdf"

DS_STATS=$TEMP_DIR"DESeq_stats_"$PROJ_ID"-"$OUTPUT_TAG".txt"

echo "There are "$COND1_NUM" samples in group A and "$COND2_NUM" samples in group B"

echo $READ_FILE_CDS $TEMP_DIR $READ_FILE_CDS
 
#DESeq1R_maker $READ_OUTPUT $KEY_OUTPUT $DS1_RESULTS $DS1_MA $OUTPUT_TAG $COND1_NUM $COND2_NUM > $DS1_RFILE 
DESeq2R_maker $READ_OUTPUT_CDS $KEY_OUTPUT $DS2_RESULTS_CDS $DS2_MA_CDS $OUTPUT_TAG $test > $DS2_RFILE_CDS 
DESeq2R_maker $READ_OUTPUT_ALL $KEY_OUTPUT $DS2_RESULTS_ALL $DS2_MA_ALL $OUTPUT_TAG $test > $DS2_RFILE_ALL 

# echo "Running DESeq1..."
# cat $DS1_RFILE | R --vanilla > $DS1_OUTPUT
echo "Running DESeq2 for CDS to generate $DS2_OUTPUT_CDS..."
cat $DS2_RFILE_CDS | R --vanilla > $DS2_OUTPUT_CDS
echo "Running DESeq2 for All Features..."
cat $DS2_RFILE_ALL | R --vanilla > $DS2_OUTPUT_ALL

echo "For DESeq2..." > $DS_STATS
up=`cat $DS2_RESULTS_CDS | awk '{if(2^$3 >= '$min_fold' && $NF <= '$max_P') print $0}' | wc -l`
dn=`cat $DS2_RESULTS_CDS | awk '{if(2^$3 <= 1/'$min_fold' && $NF <= '$max_P') print $0}' | wc -l`
echo CDS Up: $up >> $DS_STATS
echo CDS Dn: $dn >> $DS_STATS

up=`cat $DS2_RESULTS_ALL | awk '{if(2^$3 >= '$min_fold' && $NF <= '$max_P') print $0}' | wc -l`
dn=`cat $DS2_RESULTS_ALL | awk '{if(2^$3 <= 1/'$min_fold' && $NF <= '$max_P') print $0}' | wc -l`
echo ALL Up: $up >> $DS_STATS
echo ALL Dn: $dn >> $DS_STATS

echo "cp $DS2_RESULTS_CDS $DS2_RESULTS_ALL $DS2_MA_CDS $DS2_MA_ALL $DS_STATS $OUT_DIR"
cp $DS2_RESULTS_CDS $DS2_RESULTS_ALL $DS2_MA_CDS $DS2_MA_ALL $DS_STATS $OUT_DIR
#ls -lrt $OUT_DIR*
echo $OUT_DIR
