#!/bin/sh

MOC_DB="/broad/IDP-Dx_storage/MOC/MOC_DB/"

DB_FILE=`ls -lrt $MOC_DB* | tail -1 | awk '{print $9}'`

#sh /broad/IDP-Dx_storage/MOC/scripts/GS_import.sh

########################   FIELD_HEADER  ############################
FIELD_HEADER ()
{

	FILE=$1
	Q_HEAD=$2

	cat $FILE | awk '{if($1 !~ /#/) print $0}' | head -1 | sed 's/ /_/g' | awk -F"\t" -v Q_HEAD=$Q_HEAD '{	
																for(i=1; i < NF+1; i++)															
																	if($i==Q_HEAD)
																		print i
															}' 
}

#######################   FIND_FIELD  ##############################
FIND_FIELD ()
{
		
	# Look in FILE, find line in which the value in Q_COL equals Q_VAL and return the values in R_COL 
	
	Q_VAL=$1
	shift
	Q_HEAD=$1
	shift
	R_HEAD=$1
	shift
	
	
	for FILE in $@
	do

		Q_COL=`FIELD_HEADER $FILE $Q_HEAD` 
		R_COL=`FIELD_HEADER $FILE $R_HEAD` 
		
		cat $FILE | grep -v "#" | sed 's/ /_/g' | awk -F"\t" -v Q_COL=$Q_COL -v R_COL=$R_COL -v Q_VAL=$Q_VAL '{	

																			if($Q_COL==Q_VAL && R_COL != "")
																				print $R_COL
															}'  
			
	done		
														

} 
##############################


LC_F=`FIELD_HEADER $DB_FILE Library_Construction`
STATUS_F=`FIELD_HEADER $DB_FILE Status`
STATUS_F=`FIELD_HEADER $DB_FILE Status`
SAM_NUM_F=`FIELD_HEADER $DB_FILE Number_of_Samples`


ALL_LC="RNAtag-Seq 16S"

for LC in $ALL_LC
do
	echo "Process: "$LC
	awk -F"\t" -v LC_F=$LC_F -v LC=$LC -v SAM_NUM_F=$SAM_NUM_F -v total_sam=0 -v total_proj=0 '{
																										if($LC_F==LC) 
																										{
																										
																											total_sam=total_sam+$SAM_NUM_F
																											total_proj++
																										}
																										
																										print "total_samples = "total_sam
																										print "total_projects = "total_proj
																								}' $DB_FILE | tail -2
																								
	awk -F"\t" -v LC_F=$LC_F -v LC=$LC -v STATUS_F=$STATUS_F '{if($LC_F==LC)print $STATUS_F}' $DB_FILE | sort | uniq -c | sort


done