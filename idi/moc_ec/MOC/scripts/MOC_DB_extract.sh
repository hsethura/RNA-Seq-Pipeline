#!/bin/sh

DB_FILE=$1
SUBJ_HEADER=$2
SUBJ_VALUE=$3
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

DB_VALUE_FIND ()
{	
	
	
	FILE=$1
	SUBJ_H=$2
	SUBJ_V=$3
	QUERY_H=$4
	
	SUBJ_F=`FIELD_HEADER $DB_FILE $SUBJ_H`
	QUERY_F=`FIELD_HEADER $DB_FILE $QUERY_H`

	cat $FILE | grep -v "#" | awk -F"\t" -v SUBJ_F=$SUBJ_F -v QUERY_F=$QUERY_F -v SUBJ_V=$SUBJ_V -v QUERY_H=$QUERY_H '{	

														if($SUBJ_F==SUBJ_V)
															printf "%s\n", $QUERY_F
													}' 
}
 

DB_VALUE_FIND $DB_FILE $SUBJ_HEADER $SUBJ_VALUE $QUERY_HEADER
