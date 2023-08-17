#!/bin/sh



SUBJECT_DB_FILE=$1
QUERY_DB_FILE=$2
SUBJ_HEADER=$3
QUERY_HEADER=$4

##############  FUNCTIONS  #################


FIELD_HEADER ()
{
	FILE=$1
	HEADER_NAME=$2

	cat $FILE | grep -v "#" | head -1 | awk -F"\t" -v HEADER_NAME=$HEADER_NAME '{	
																for(i=1; i < NF+1; i++)															
																	if($i==HEADER_NAME)
																		print i
															}' 
}


 
	
SUBJ_F=`FIELD_HEADER $SUBJECT_DB_FILE $SUBJ_HEADER`
QUERY_F=`FIELD_HEADER $QUERY_DB_FILE $QUERY_HEADER`
SUBJ_REF_F=`FIELD_HEADER $SUBJECT_DB_FILE Bacterial_reference`

typeset -i NL
typeset -i i

NL=`cat $SUBJECT_DB_FILE | grep -v "#" | wc -l`
p=1

cat $SUBJECT_DB_FILE | grep -v "#" | head -1 | awk -F"\t" -v p=$p -v SUBJ_F=$SUBJ_F '{ 				
																									for(i=1; i < NF+1; i++) 
																										printf "%s\t", $i
																						}'
cat $QUERY_DB_FILE | grep -v "#" | head -1 | awk -F"\t" -v p=$p -v SUBJ_F=$SUBJ_F '{ 				
							for(i=1; i < NF+1; i++) 
								printf "%s\t", $i
				}'
				
printf "\n"																			

while [ $p -le $NL ]
do
	
	SUBJ_V=`cat $SUBJECT_DB_FILE | grep -v "#" | sed 1d | awk -F"\t" -v p=$p -v SUBJ_F=$SUBJ_F '{
																					if(NR==p)
																						print $SUBJ_F
																				}'`
	SUBJ_REF=`cat $SUBJECT_DB_FILE | grep -v "#" | sed 1d | awk -F"\t" -v p=$p -v SUBJ_REF_F=$SUBJ_REF_F '{
																					if(NR==p)
																						print $SUBJ_REF_F
																				}' | head -1`

	cat $SUBJECT_DB_FILE | grep -v "#" | sed 1d | awk -F"\t" -v p=$p '{
				
				if(NR==p)
				{
					for(i=1; i < NF+1; i++) 
						printf "%s\t", $i
				}
			}' 

 	cat $QUERY_DB_FILE | grep -v "#" | sed 's/_'$SUBJ_REF'//' | sed 1d | awk -F"\t" -v p=$p -v SUBJ_V=$SUBJ_V -v QUERY_F=$QUERY_F '{

																						if($QUERY_F==SUBJ_V)
																							for(i=1; i < NF+1; i++) 
																								printf "%s\t", $i
																					}'
	printf "\n"																			
	
	p=`expr $p + 1`
done


