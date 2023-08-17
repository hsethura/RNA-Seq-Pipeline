#!/bin/sh


MOC_ID=$1
QUERY_DB_PATH=$2
QUERY_HEADER=$3

ID=`id | cut -d"(" -f2 | cut -d")" -f1`
KEY_DIR="/broad/IDP-Dx_storage/MOC/Key_files/"
JOIN_DIR="/broad/IDP-Dx_storage/MOC/MKJoin_files/"
KEY_FILE=$KEY_DIR"/"$MOC_ID"_key.txt"
JOIN_FILE=$JOIN_DIR"/"$MOC_ID"_KeyMetrics.txt"
KEY_HEADER="Sample_ID"

mkdir -p $JOIN_DIR
echo "" | sed 1d > $JOIN_FILE

echo $JOIN_FILE


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


 
KEY_F=`FIELD_HEADER $KEY_FILE $KEY_HEADER`
KEY_REF_F=`FIELD_HEADER $KEY_FILE Bacterial_reference`
KEY_PID_F=`FIELD_HEADER $KEY_FILE Project_ID`
PATH_F=`FIELD_HEADER $KEY_FILE Path_to_SeqFile`

typeset -i NL
typeset -i i

NL=`cat $KEY_FILE | grep -v "#" | wc -l`
p=1


while [ $p -le $NL ]
do
	
	KEY_SAMPID=`cat $KEY_FILE | grep -v "#" | sed 1d | awk -F"\t" -v p=$p  '{
																					if(NR==p)
																					{
																						if($'$KEY_F' != "")
																							print $'$KEY_F'
																						if($'$KEY_F' == "")
																							print ""
																					}
																				}'`
	KEY_REF=`cat $KEY_FILE | grep -v "#" | sed 1d | awk -F"\t" -v p=$p '{
																					if(NR==p)
																					if($'$KEY_REF_F' ~ /;/)
																					{
																						split($'$KEY_REF_F', ar, ";")
																						print ar[2]
																					}
																					else
																						print $'$KEY_REF_F'
																				}' | head -1`

	KEY_PROJID=`cat $KEY_FILE | grep -v "#" | sed 1d | awk -F"\t" -v p=$p '{
																					if(NR==p)
																						print $'$KEY_PID_F'
																				}' | head -1`

	if [ $KEY_SAMPID != "" ];then
		#QUERY_DB_FILE=$QUERY_DB_PATH"/"$ID"/"$KEY_PROJID"/"$KEY_PROJID"_metrics.txt"
		QUERY_DB_FILE=$QUERY_DB_PATH"/"$KEY_PROJID"/"$KEY_PROJID"_metrics.txt"

		MET_REP_F=`FIELD_HEADER $QUERY_DB_FILE "For_replicon..."`

		ALL=`cat $QUERY_DB_FILE | awk -F"\t" '{ print $'$MET_REP_F'}' | grep "ALL" | wc -l`
	
		if [ $p == 1 ];then
			
			cat $KEY_FILE | grep -v "#" | head -1 | awk -F"\t" -v p=$p -v KEY_F=$KEY_F '{ 				
																							for(i=1; i < '$PATH_F'+1; i++) 
																								printf "%s\t", $i

																				}' >> $JOIN_FILE
			
			QUERY_F=`FIELD_HEADER $QUERY_DB_FILE $QUERY_HEADER`
			cat $QUERY_DB_FILE | grep -v "#" | head -1 | awk -F"\t" -v p=$p -v KEY_F=$KEY_F '{ 				
								for(i=1; i < NF+1; i++) 
									printf "%s\t", $i
					}' >> $JOIN_FILE
			
	
			printf "\n"	>> $JOIN_FILE																		
		fi
		cat $KEY_FILE | grep -v "#" | sed 1d | awk -F"\t" -v p=$p '{
			
				if(NR==p)
				{
					for(i=1; i < '$PATH_F'+1; i++) 
						printf "%s\t", $i
				}
			}' >> $JOIN_FILE

		cat $QUERY_DB_FILE | grep -v "#" | sed 's/_'$KEY_REF'//' | sed 1d | awk -F"\t" -v p=$p -v KEY_SAMPID=$KEY_SAMPID -v QUERY_F=$QUERY_F -v ALL=$ALL '{

																							
																								if($QUERY_F==KEY_SAMPID && (ALL=="0" || $'$MET_REP_F'=="ALL"))
																									for(i=1; i < NF+1; i++) 
																										printf "%s\t", $i
																							}' >> $JOIN_FILE
		printf "\n"	>> $JOIN_FILE																	
	fi
	p=`expr $p + 1`
done


