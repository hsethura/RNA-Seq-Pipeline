#!/bin/sh

CONFIG_FILE=$1
shift
KEY_ID=$1
shift
PROJ_ID=$1
shift
ACC=$1
shift
COND1_ID=$1
shift
COND2_ID=$1
shift
COND1_TAG=$1
shift
COND2_TAG=$1
shift

source "/broad/software/scripts/useuse"
reuse -q UGER
reuse -q Java-1.8
reuse -q R-3.2

ID=`id | cut -d"(" -f2 | cut -d")" -f1`


##############  FUNCTIONS  #################

# Pulling out paths, header names from config file 
	config_read ()
	{
	
		CONFIG_FILE=$1
		VAR_NAME=$2
		cat $CONFIG_FILE | grep $VAR_NAME":" | awk '{print $2}'
	
	
	
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
	
	echo "library(DESeq2)"
	echo "counts = read.table(\""$READ_OUTPUT"\", header=T, row.names = 1, sep = \"\t\", quote = \"\", comment.char = \"\")" 
	echo "key = read.table(\""$KEY_OUTPUT"\", header=T)" 
	echo "rownames(key) = key[,1]"
	echo "des = DESeqDataSetFromMatrix(counts, key, ~Group)" 
	echo "des = DESeq(des, modelMatrixType = \"expanded\")" 
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
	DS2_MA=$4
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



KEY_PATH=`config_read $CONFIG_FILE Key_base`
RESULTS_PATH=`config_read $CONFIG_FILE Results_path`
DESEQ_PATH=`config_read $CONFIG_FILE DESeq_path`
CGID_NAME=`config_read $CONFIG_FILE CGID_NAME`

KEY_FILE=$KEY_PATH"/"$KEY_ID"_key.txt"
TEMP_KEY_FILE=$KEY_PATH"/"$PROJ_ID"_temp_key.txt"


TAG=$COND1_TAG"_vs_"$COND2_TAG
PROJ_DIR=$RESULTS_PATH"/"$PROJ_ID"/"
DS_DIR=$DESEQ_PATH"/"$ID"/"$PROJ_ID"/DESeq/"
mkdir -p $DS_DIR
READ_FILE_TEMP=$PROJ_DIR"/"$PROJ_ID"_"$ACC".tsv"
READ_FILE=$DS_DIR"/"$PROJ_ID"=REP_Genes_"$ACC"_Reads.txt"
DS1_RFILE=$DS_DIR$PROJ_ID"-"$TAG"_DS1.R"
DS2_RFILE=$DS_DIR$PROJ_ID"-"$TAG"_DS2.R"
TEMP_OUTFILE1=$DS_DIR"/temp1.txt"
TEMP_OUTFILE2=$DS_DIR"/temp2.txt"
TEMP_OUTFILE3=$DS_DIR"/temp3.txt"


echo "cat $READ_FILE_TEMP | grep -v tRNA: | grep -v IGR: | grep -v AS_ | grep -v rRNA: > $READ_FILE"

if [ -s $READ_FILE_TEMP ];then
	cat $READ_FILE_TEMP | grep -v tRNA: | grep -v IGR: | grep -v AS_ | grep -v rRNA: > $READ_FILE
else
	echo $READ_FILE_TEMP" not found"
fi

ls -lrt $KEY_FILE
ls -lrt $READ_FILE


if [ ! -s $READ_FILE ];then

	echo $READ_FILE" not found!"

fi
if [ ! -s $KEY_FILE ];then

	echo $KEY_FILE" not found!"

fi


min_fold=`extract_option -f 2 1 $@`
max_P=`extract_option -P 0.05 1 $@`
include=`extract_option -i - 100 $@`
compare=`extract_option -c - 2 $@`
omit=`extract_option -o ABCDEFGHIJKLMNOP 100 $@`

			
# IDF=`cat $KEY_FILE | head -3 | tail -1 | awk -F"\t" '{
# 														for(i=1; i < NF+1; i++)
# 															if($i=="Sample_ID")
# 																print i
# }'`
# 
# C1=`cat $KEY_FILE | awk -v m=$m -v IDF=$IDF -v COND1_ID=$COND1_ID '{
# 																			y=split(COND1_ID, ar, ",")
# 																			for(i=1; i < y+1; i++)
# 																				if(NR==(ar[i]+3))
# 																					print $IDF
# 																	
# 																	}'`
# C2=`cat $KEY_FILE | awk -v m=$m -v IDF=$IDF -v COND2_ID=$COND2_ID '{
# 																			y=split(COND2_ID, ar, ",")
# 																			for(i=1; i < y+1; i++)
# 																				if(NR==(ar[i]+3))
# 																					print $IDF
# 																	
# 																	}'`

ALL_COND=`echo $COND1_ID $COND2_ID | sed 's/,/ /g'`

echo $ALL_COND

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

KEY_OUTPUT=$DS_DIR$PROJ_ID"-"$TAG"_key.txt"
READ_OUTPUT=$DS_DIR$PROJ_ID"-"$TAG"_reads.txt"

cat $TEMP_OUTFILE1 > $READ_OUTPUT

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


tag2=$tag"_"$numA"vs"$numB
DS1_RESULTS=$DS_DIR"DESeq1_results_"$PROJ_ID"-"$OUTPUT_TAG".txt"
DS2_RESULTS=$DS_DIR"DESeq2_results_"$PROJ_ID"-"$OUTPUT_TAG".txt"
DS1_OUTPUT=$DS_DIR"DESeq1_out"$PROJ_ID"-"$OUTPUT_TAG".txt"
DS2_OUTPUT=$DS_DIR"DESeq2_out_"$PROJ_ID"-"$OUTPUT_TAG".txt"
DS1_MA=$DS_DIR"DESeq1_"$PROJ_ID"-"$OUTPUT_TAG"_MA.pdf"
DS2_MA=$DS_DIR"DESeq2_"$PROJ_ID"-"$OUTPUT_TAG"_MA.pdf"

echo "There are "$COND1_NUM" samples in group A and "$COND2_NUM" samples in group B"


DESeq1R_maker $READ_OUTPUT $KEY_OUTPUT $DS1_RESULTS $DS1_MA $OUTPUT_TAG $COND1_NUM $COND2_NUM > $DS1_RFILE 
DESeq2R_maker $READ_OUTPUT $KEY_OUTPUT $DS2_RESULTS $DS2_MA $OUTPUT_TAG > $DS2_RFILE 

echo "Running DESeq1..."
cat $DS1_RFILE | R --vanilla > $DS1_OUTPUT
echo "Running DESeq2..."
cat $DS2_RFILE | R --vanilla > $DS2_OUTPUT


ls -lrt $DS1_RESULTS
ls -lrt $DS2_RESULTS
ls -lrt $DS2_RFILE
ls -lrt $DS1_MA
ls -lrt $DS2_MA

echo "For DESeq1..."
up=`cat $DS1_RESULTS | awk '{if(2^$6 >= '$min_fold' && $8 <= '$max_P') print $0}' | wc -l`
dn=`cat $DS1_RESULTS | awk '{if(2^$6 <= 1/'$min_fold' && $8 <= '$max_P') print $0}' | wc -l`
echo Up: $up
echo Dn: $dn

echo "For DESeq2..."
up=`cat $DS2_RESULTS | awk '{if(2^$3 >= '$min_fold' && $NF <= '$max_P') print $0}' | wc -l`
dn=`cat $DS2_RESULTS | awk '{if(2^$3 <= 1/'$min_fold' && $NF <= '$max_P') print $0}' | wc -l`
echo Up: $up
echo Dn: $dn

ls -lrt $DS_DIR*