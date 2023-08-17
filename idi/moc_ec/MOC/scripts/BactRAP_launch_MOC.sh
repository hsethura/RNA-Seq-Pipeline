#!/bin/sh

source "/broad/software/scripts/useuse"
reuse -q UGER
reuse -q Java-1.6


EXP_ID=$1
shift 
options=$@
 
options=$options" -"

OUT_DIR="/broad/IDP-Dx_storage/Bact_RNASeq_results/"
SCRIPT_PATH="/broad/IDP-Dx_storage/BactRAP/scripts/"
ACC_PATH="/broad/IDP-Dx_storage/NCBI_files2/"
INPUT_DIR="/broad/IDP-Dx_storage/MOC/key_files/"
DB_FILE="/broad/IDP-Dx_storage/BactRAP/scripts/all_NCBI_DB.txt"
ANN_NUC_PATH="/broad/IDP-Dx_storage/BactRAP/genome_annotations/"

INPUT_FILE=$INPUT_DIR$EXP_ID"_key.txt"
TEMP_DIR="/broad/hptmp/Dx/"$EXP_ID"/"
RUN_MERGEDIN_FILE=$TEMP_DIR$EXP_ID"_input_temp.txt"
ALIGN_COMMAND_FILE=$TEMP_DIR$EXP_ID"_align_commands.txt"

echo $INPUT_FILE

cat $INPUT_FILE | sed 's/ /_/g' > temp_input.txt
INPUT_FILE="temp_input.txt"

##############################################  START OF FUNCTIONS  ######################################
########   function for testing UGER status
JOB_test ()
{	
	
	running=`qstat | sed 1d | sed 1d | awk '{print $1}'`
	echo "These jobs are running: " $running
	echo "Waiting for the following jobs to finish: " $@

	while :
	do 
		running=`qstat | sed 1d | sed 1d | awk '{print $1}'`
		jobs=`echo $running $@ | awk '{
									for(i=1; i < NF+1; i++)
										printf "%s\n", $i
									}' | sort | uniq -c | awk '{if($1 == 2) print $0}' | wc -l`
		if [ $jobs -eq 0 ];then
			break
		fi
	done 
}

#########  function for extracting options from command line
 
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
			{
				print $(i+1)
			}
		
		}

	}' | tail -1 
}

########## function for extracting field names from key headers

find_field()
{
	INPUT_FILE=$1
	field_title=$2
	req_field=$3
	
	cat $INPUT_FILE | head -3 | tail -1 | awk -F"\t" -v f=$field_title -v req=$req_field '{
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

##############################################  END OF FUNCTIONS  ######################################


############  Extract field IDs from key headers  #################

flow_field=`find_field $INPUT_FILE Flowcell N`
lane_field=`find_field $INPUT_FILE Lane N`
dir_field=`find_field $INPUT_FILE Directory N`
ID_field=`find_field $INPUT_FILE Sample_ID Y`
Acc_field=`find_field $INPUT_FILE RefSeq_accession Y`
PE=`find_field $INPUT_FILE Read_Pairing Y`
Condition_1=`find_field $INPUT_FILE Condition_1 N`
Resistance=`find_field $INPUT_FILE Resistance? N`
Inline_Barcode_seq=`find_field $INPUT_FILE Inline_Barcode_seq N`
seq=`find_field $INPUT_FILE Sequencer N`
P7=`find_field $INPUT_FILE P7_Index N`
LC=`find_field $INPUT_FILE LC_strategy N`
Project_ID=`find_field $INPUT_FILE Project_ID N`

echo $flow_field flow_field
echo $lane_field lane_field
echo $dir_field dir_field
echo $ID_field ID_field
echo $bcParsed_files bcParsed_files
echo $Acc_field Acc_field
echo $PE PE
echo $RF_field RF_field
echo $Genus_field Genus_field
echo $Strain Strain

f=$ID" "$Acc_field" "$PE" "$seq" "$Project_ID

######################  test if required field is missing ########################
for m in $f
do
	exit=`echo $m | cut -d";" -f2`
	title=`echo $m | cut -d";" -f1`
	
	if [ $exit == "exit" ];then
		
		echo "Required field "$title" is missing!"
		exit
	
	fi
done

############################  extract options  ##############################

UGER=`extract_option -UGER Y $options`
SCRIPT_PATH=`extract_option -SCRIPT_PATH $SCRIPT_PATH $options`
OUT_DIR=`extract_option -OUT_DIR $OUT_DIR $options`
GL_ow=`extract_option -GL_ow N $options`
AH_ow=`extract_option -AH_ow - $options`
q=`extract_option -q short $options`
ANN_OW=`extract_option -ANN_OW N $options`
MAIL=`extract_option -MAIL Y $options`
UGER_HEAD=`extract_option -UGER_HEAD 100000000 $options`
UGER_TAIL=`extract_option -UGER_TAIL 100000000 $options`
COMP=`extract_option -COMP N $options`
ADD_GFF=`extract_option -ADD_GFF N $options`
BEST_ACC=`extract_option -BEST_ACC N $options`
EXPID=`extract_option -EXPID $EXP_ID $options`
ACC_PATH=`extract_option -ACC_PATH $ACC_PATH $options`

RUN_BCP=`extract_option -RUN_BCP N $options`
RUN_BACT=`extract_option -RUN_BACT N $options`

#########################  Set directory and file paths  #######################

TEMP_DIR="/broad/hptmp/Dx/"$EXPID"/"
EXP_DIR=$OUT_DIR$EXPID"/"
LAUNCH_FILE=$EXP_DIR"_LAUNCH_FILE.txt"
COMMAND_FILE=$EXP_DIR"commands.txt"
ANALYSIS_FILE=$EXP_DIR"analysis_commands.txt"
DESEQ2_KEY=$EXP_DIR$EXPID"_DES2_key.txt"
UGER_COMMAND=$EXP_DIR"UGER.txt"
JOBS_FILE=$EXP_DIR"UGER_jobs.txt"

mkdir $EXP_DIR 
mkdir $TEMP_DIR 

SEQ=`cat $INPUT_FILE | sed 1d | awk -F"\t" -v seq=$seq '{print $seq}' | uniq`

echo "" | sed 1d > $JOBS_FILE

########################## Generate command lines for BactRAP #########################

echo "Generating command lines for BactRAP..."


echo "Sample_ID	Strain_Abx_tp" > temp_key.txt

# pull out the information from the bcParsed_files INPUT_FILE to set up command lines for pipeline


cat $INPUT_FILE | awk -F"\t" -v EXPID=$EXPID -v EXP_DIR=$EXP_DIR -v LSF=$LSF -v TEMP_DIR=$TEMP_DIR -v SCRIPT_PATH=$SCRIPT_PATH -v OUT_DIR=$OUT_DIR -v options="$options" \
								 -v merged_file=$merged_file -v ID_field=$ID_field -v Acc_field=$Acc_field -v PE=$PE -v Inline_Barcode_seq=$Inline_Barcode_seq \
								 -v Condition_1=$Condition_1 -v Resistance=$Resistance -v Time=$Time -v LC=$LC -v Inline_Barcode_seq=$Inline_Barcode_seq -v Project_ID=$Project_ID '{

						if(NR > 3 && $1 != "")
						{	
							
							bc=substr($Inline_Barcode_seq, 1, 8)
							merged_file=TEMP_DIR"/merged_"bc"_R1.fastq"
							
							ANALYSIS_ID="ID_"$ID_field"=PE_"$PE"End=" 
							OUT_FILE=ANALYSIS_ID"_out.txt"
						
							if($PE=="Paired")
								PE_ID="Y"
							if($PE=="Single")
								PE_ID="N"
							if($PE!="Single" && $PE!="Paired")
								PE_ID=$PE
							
							if($Acc_field ~/;/)
							{	
								split($Acc_field, ar, ";")
								ACC=ar[2]
							}
							else
								ACC=$Acc_field
								
							
							print "sh "SCRIPT_PATH"BactRAP.sh "merged_file" "$Project_ID" "$Acc_field" -ANALYSIS_ID "ANALYSIS_ID" -SAMPLE_ID "$ID_field" -SCRIPT_PATH "SCRIPT_PATH" -OUT_DIR "OUT_DIR" -PE "PE_ID" -LC TS "options > "temp_commands.txt"
							print "qsub -l h_vmem=6g -b y sh "SCRIPT_PATH"BactRAP.sh "merged_file" "$Project_ID" "$Acc_field" -ANALYSIS_ID "ANALYSIS_ID" -SAMPLE_ID "$ID_field" -SCRIPT_PATH "SCRIPT_PATH" -OUT_DIR "OUT_DIR" -PE "PE_ID" -LC TS "options > "temp_UGER.txt"
						} 
		
					}'  
 
  
echo "temp_commands.txt"
echo "temp_UGER.txt"