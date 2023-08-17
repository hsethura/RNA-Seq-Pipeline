#!/bin/sh


MOC_ID=$1
shift
CONFIG_FILE=$1
shift
USERID=$1
shift

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"


ID=`id | cut -d"(" -f2 | cut -d")" -f1`
KEY_DIR="/broad/IDP-Dx_storage/MOC/Key_files/"
KEY_FILE=$KEY_DIR"/"$MOC_ID"_key.txt"
JL_SCRIPTS_DIR="/broad/IDP-Dx_storage/MOC/scripts/"
NB_SCRIPTS_DIR="/broad/IDP-Dx_work/nirmalya/research/edgeR_scripts/"


if [ $USERID == "N" ];then
	ID=""
fi

echo $ID


TEMP_FILE=~/"temp_"$MOC_ID"_DSlaunch.txt"

echo $KEY_FILE

cat $KEY_FILE | sed 's/, /,/g' | sed 's/ /_/g' | sed 's*/*_*g' > $TEMP_FILE



# typeset -i SAMP_ID_FIELD
# typeset -i CG_ID_FIELD


##############  FUNCTIONS  #################

# Pulling out paths, header names from config file 
config_read ()
{

	CONFIG_FILE=$1
	VAR_NAME=$2
	cat $CONFIG_FILE | grep $VAR_NAME":" | awk '{print $2}'



}

FIELD_HEADER ()
{

	KEY_FILE=$1
	HEADER_NAME=$2

	cat $KEY_FILE | grep -v "#" | head -1 | awk -F"\t" -v HEADER_NAME=$HEADER_NAME '{	
																if(NR==1)
																	print "0"

																for(i=1; i < NF+1; i++)															
																	if($i==HEADER_NAME)
																		print i
															}' | tail -1 
}

FIND_FIELD ()
{
	FILE=$1
	shift
	GCID=$1
	shift
	SIF=$1
	shift
	
	
	for CGF in $@
	do
	
		cat $FILE | grep -v "#" | awk -F"\t" -v GCID=$GCID -v CGF=$CGF -v SIF=$SIF -v EXCLUDE=$EXCLUDE '{	

														
															y=split(GCID, ar, ";")
															
															for(i=1; i < y+1; i++)
															{
																if($CGF == ar[i])
																{
																	if($SIF ~ /;/)
																	{	
																		split($SIF, ar, ";")
																		printf ",%s,", ar[2]
																	}
																	if($SIF !~ /;/)
																		printf "%s,", $SIF
																}	
															}
														}'
	done 

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

EXCLUDE=`extract_option -EXC - 100 $@`
INCLUDE=`extract_option -INC - 100 $@`

EDGR_PATH=`config_read $CONFIG_FILE DE_path`
RESULTS_PATH=`config_read $CONFIG_FILE Results_path`

SAMP_ID_FIELD=`FIELD_HEADER $KEY_FILE Sample_ID`
CG_ID_FIELD=`FIELD_HEADER $KEY_FILE CG_ID`
PID_FIELD=`FIELD_HEADER $KEY_FILE Project_ID`
ACC_FIELD=`FIELD_HEADER $KEY_FILE Bacterial_reference`
PAIRS_FIELD=`FIELD_HEADER $KEY_FILE CG_Pairs`

ALL_PAIRS=`cat $TEMP_FILE | awk -F"\t" '{if(NR>5) print $0}' | awk -F"\t" -v PAIRS_FIELD=$PAIRS_FIELD '{print $PAIRS_FIELD}' | sed '/^$/d' | sed '/^$/d' | sort | uniq`

# echo $PAIRS_FIELD
# echo $ALL_PAIRS
if [ $PAIRS_FIELD == "0" ];then

	echo "No CG_Pairs field found..."
	exit
fi

echo $ALL_PAIRS 
echo $ALL_PAIRS | wc -w

for m in $ALL_PAIRS
do
	
	GCIDA=`echo $m | cut -d"," -f1 |sed 's/:/,/g'`
	GCIDB=`echo $m | cut -d"," -f2 | sed 's/:/,/g'`

	A=`FIND_FIELD $TEMP_FILE $GCIDA $SAMP_ID_FIELD $CG_ID_FIELD $SAMP_ID_FIELD | sed 's/,$//g' | tr ',' '\n' | sed '/^$/d' | sort -u | tr '\n' ','` 
	B=`FIND_FIELD $TEMP_FILE $GCIDB $SAMP_ID_FIELD $CG_ID_FIELD $SAMP_ID_FIELD | sed 's/,$//g' | tr ',' '\n' | sed '/^$/d' | sort -u | tr '\n' ','` 
	PROJ_ID=`FIND_FIELD $TEMP_FILE $GCIDA $PID_FIELD $CG_ID_FIELD $SAMP_ID_FIELD | tr ',' '\n' | sort -u` 
	ACC=`FIND_FIELD $TEMP_FILE $GCIDA $ACC_FIELD $CG_ID_FIELD $CG_ID_FIELD $SAMP_ID_FIELD | sed 's/,,/,/g' | tr ',' '\n' | sed '/^$/d'  | sort -u` 

	PROJ_DIR=$RESULTS_PATH"/"$ID"/"$PROJ_ID"/"
	EDGR_DIR=$EDGR_PATH"/"$ID"/"$PROJ_ID"/DE/EdgeR/"
	READ_FILE=$PROJ_DIR"/"$PROJ_ID"_"$ACC"_counts.tsv"
	GENE_READ_FILE=$PROJ_DIR"/"$PROJ_ID"_"$ACC"_genes_counts.tsv"
	EDGER_READ_FILE=$PROJ_DIR"/"$PROJ_ID"_"$ACC"_edgeR_counts.tsv"

	#cat $READ_FILE | awk '{if($1 ~ /^CDS:/ || $1 ~ /^ncRNA:/ || NR==1) print $0}' > $GENE_READ_FILE
	cat $READ_FILE | tr -d \'\" > $EDGER_READ_FILE
	
	echo $EDGER_READ_FILE
		
	mkdir $PROJ_DIR
 
 	echo "#########################"
	echo $PID_FIELD
	echo "A" $A 
	echo "B" $B 
	echo "GCIDA  "$GCIDA
	echo "GCIDB  "$GCIDB
	echo "ACC  "$ACC
	echo "READ_FILE	"$READ_FILE	
	echo "EDGER_READ_FILE "$EDGER_READ_FILE
	echo "EDGR_DIR "$EDGR_DIR
	echo "#########################"
		

	if [ -s $READ_FILE ];then
		
		echo "sh $JL_SCRIPTS_DIR"DESeq_launch5.sh" $CONFIG_FILE $MOC_ID $PROJ_ID $ACC $A $B $GCIDA $GCIDB $USERID -test WALD"

		sh $JL_SCRIPTS_DIR"DESeq_launch5.sh" $CONFIG_FILE $MOC_ID $PROJ_ID $ACC $A $B $GCIDA $GCIDB $USERID -test WALD 

		ER_A=`echo $A | sed 's/,/;/g'`
		ER_B=`echo $B | sed 's/,/;/g'`
	
		echo $B
	
	echo "$NB_SCRIPTS_DIR"edgeR_main.R" -i $EDGER_READ_FILE -c \"$ER_A\" -n \"$ER_B\" -t lr -o $EDGR_DIR --tag $GCIDA"_"$GCIDB"_lr""
		$NB_SCRIPTS_DIR"edgeR_main.R" -i $EDGER_READ_FILE -c \"$ER_B\" -n \"$ER_A\" -t lr -o $EDGR_DIR --tag $GCIDA"_"$GCIDB"_lr"

	fi
done


