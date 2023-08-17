#!/bin/sh

REF_PATH=$1"/"
shift
OUT_DIR=$1"/"
shift
TEMP_DIR=$1"/"
shift
ADD5=$1
shift
ADD3=$1
shift
ACC=$1

mkdir $TEMP_DIR"/"
mkdir $OUT_DIR"/"

##############

fna_extract ()
{

	TEMP_FNA=$1
	REP_PARSED=$2
	GFF_REP=$3
	
	typeset -i NUM_MATCH
	NUM_MATCH=`	awk -v flag="N" -v GFF_REP=$GFF_REP -v REP_PARSED=$REP_PARSED '{
			if($0 ~ /'$REP_PARSED'/ && $1 ~ />/ )
			{	
				print ">"GFF_REP
			}
		}' $TEMP_FNA | wc -l`

	
	cat $TEMP_FNA | tr '\r' '\n' | awk -v flag="N" -v GFF_REP=$GFF_REP -v REP_PARSED=$REP_PARSED -v NUM_MATCH=$NUM_MATCH '{
	
					
					if($1 ~ />/ && flag=="Y")
						flag="N"
					if(flag=="Y" && $1 !~ />/)
 						print $1
					if ($1 ~ />/ )
					{
						if(NUM_MATCH==1)
						{
					
							if($0 ~ /'$REP_PARSED'/)
							{	
								print ">"GFF_REP
								flag="Y"
							}
						}
						if(NUM_MATCH > 1)
						{
							if($1 == ">"REP_PARSED)
							{	
								print ">"GFF_REP
								flag="Y"
							}
						}
					}
	}' | sed '/^$/d' 

}

############################

gff_info_parse ()
{
	LINE=$1
	shift
	F=$1
	shift
		
	for KEY in $@
	do
		echo $LINE | awk -v KEY=$KEY '{	
	
						y=split($'$F', ar, ";")
						for(i=0; i < y+1; i++)
						{
							split(ar[i], pr, "=")
							if(pr[1]==KEY)
								print pr[2]

						}
					}'
	done
}

# input files
GFF_FILE=$REF_PATH"/"$ACC".gff"
FNA_FILE=$REF_PATH"/"$ACC".fna"

# output files
ALL_FILE=$OUT_DIR"/"$ACC"_ALL.gff"
GENE_FILE=$OUT_DIR"/"$ACC"_GENES.gff"
PARSED_FNA=$OUT_DIR$ACC".fna"
TEMP_FNA=$TEMP_DIR"/"$ACC"_temp.fna"


cat $GFF_FILE | tr '\r' '\n' | grep "#" > $ALL_FILE
cat $FNA_FILE | sed -e 's/(//g' -e 's/)//g' -e 's/|//g' > $TEMP_FNA

echo "" | sed 1d > $GENE_FILE
echo "" | sed 1d > $PARSED_FNA 


ALL_GFF_REPS=`cat $GFF_FILE | sed 1d | sed -e 's/ /_/g' -e 's/(//g' -e 's/)//g' -e 's/|//g' | awk '{if($1 ~ /FASTA/) exit; print $1}' |  sed '/^$/d' | sort -r | uniq | grep -v "#" `
echo $ALL_GFF_REPS


for REP in $ALL_GFF_REPS
do	
	echo "REP: " $REP
	
	GFF_ACC_FILE=$TEMP_DIR"/"$REP"_rep.gff"
	GENE_ACC_FILE=$TEMP_DIR"/"$REP"_rep_GENES.gff"
	IGR_ACC_FILE=$TEMP_DIR"/"$REP"_rep_IGR.gff"
	
	cat $GFF_FILE | grep -v "#" | sed -e 's/ /_/g' -e '/^$/d' -e 's/)//g' -e 's/(//g' -e 's/|//g' -e 's/Coding_gene/CDS/g' | sort -k4n | tr '\r' '\n' | awk -F"\t" -v REP=$REP '{if($3 !~ /region/ && $3 !~ /ource/ && $1 == REP) print $0}' > $GFF_ACC_FILE
	
	echo $GFF_ACC_FILE	

	ALL_COMBOS=`cat $GFF_ACC_FILE | awk '{print $1";"$4";"$5";"$7"\t"$4}' | sort -k2n | uniq | awk '{print $1}' ` 

	echo "Parsing $GFF_ACC_FILE to generate $GENE_ACC_FILE..."
		
	typeset -i NL	
	NL=`cat $GFF_ACC_FILE | wc -l`
	
	for COMBO in $ALL_COMBOS
	do
		REP=`echo $COMBO | cut -d";" -f1`
		START=`echo $COMBO | cut -d";" -f2`
		END=`echo $COMBO | cut -d";" -f3`
		DIR=`echo $COMBO | cut -d";" -f4`

		awk -F"\t" -v REP=$REP -v START=$START -v DIR=$DIR -v TAG="" -v INFO="" -v NL=$NL -v FOUND="N" '{
										
											#print $1"|"REP"|"$4"|"'$START'"|"$5"|"'$END'"|"$7"|"DIR 
											LEN='$END'-'$START'
											if($1==REP && $4=='$START' && $5=='$END' && $7==DIR) 
											{	
												F_START=$4
												F_END=$5
												F_DIR=$7

												FOUND="Y"
												if(TAG=="")
													TAG=$3
												else
													TAG=TAG","$3
										
												if(INFO=="")
													INFO=$9
												else
													INFO=INFO","$9
											}
									
											if(($4 != '$START' &&  FOUND=="Y") || NR==NL) 
											{	
												y=split(TAG, ar, ",")
												if(y > 1)
												{	
													TEMP_TAG=""
													for(i=0; i < y+1; i++)
													{
														if(ar[i] == "CDS")
															TEMP_TAG="CDS"
														if(ar[i] == "tRNA")
															TEMP_TAG="tRNA"
														if(ar[i] == "rRNA")
															TEMP_TAG="rRNA"
														if(ar[i] == "miscRNA" || ar[i] == "misc_RNA" || ar[i] == "ncRNA" || ar[i] == "tmRNA")
															TEMP_TAG="ncRNA"
													}
													if(TEMP_TAG=="")
														TAG=ar[y]
													else
														TAG=TEMP_TAG
												}					
												print INFO 
												y=split(INFO, ar, ";")
												LC="-"
												GENE="-"
												PROD="-"
												ID="-"
												for(i=1; i < y+1; i++)
												{
													split(ar[i], pr, "=")
													if(pr[1]=="locus_tag")	
														LC=pr[2] 
													if(pr[1]=="gene" || pr[1]=="Name" )
														GENE=pr[2] 
													if(pr[1]=="product")
														PROD=pr[2] 
													if(pr[1]=="ID")
														ID=pr[2] 
												}
												
												if(TAG=="CDS")
												{
													if($7=="+")
													{
														F_START=$4-'$ADD5'
														F_END=$5+'$ADD3'
													}
													if($7=="-")
													{
														F_START=$4-'$ADD3'
														F_END=$5+'$ADD5'
													}
													if(F_START < 0)
														F_START=1	
												}
												
												
												FULL_INFO="full_tag="TAG":"LC":"GENE":"PROD":"ID":"DIR":"LEN
												
												print REP"	"$2"	"TAG"	"F_START"	"F_END"	.	"F_DIR"	"$8"	"FULL_INFO
												exit
											}

										}' $GFF_ACC_FILE

	done | sed 's/full_tag=transcript/full_tag=ncRNA/g' | sed 's/full_tag=tmRNA/full_tag=ncRNA/g' | sed 's/type=tmRNA/type=ncRNA/g' 
	#> $GENE_ACC_FILE 
	
	
	echo "Parsing $GENE_ACC_FILE to generate $IGR_ACC_FILE..."

	echo "" | sed 1d > $IGR_ACC_FILE

	cat $GENE_ACC_FILE | awk -F"\t" -v last_end=1 -v last_name=="ORIGIN" -v last_tag="ORIGIN" -v last_dir="|" '{

								$4=$4-1;$5=$5+1
								if($4 < 0)
									$4=0								
							
							
								y=split($9, ar, ";")
								for(i=0; i < y+1; i++)
								{
									split(ar[i], pr, "=")
									if(pr[1]=="type")
									{	
										type_name=pr[2]
										type_field=ar[i]
									}
									if(pr[1]=="full_tag")
									{
										tag_name=pr[2]
										tag_field=ar[i]
									}
								}

								dir=$7
								strand="+"
							
								if(NR >= 1)
								{
							 
									$3="IGR"
									if(last_end < $4)
									{	
										gsub(":", ",", last_tag)
										gsub(":", ",", tag_name)

										print $1, $2, $3, last_end, $4, $6, strand, $8, "type=IGR;full_tag=IGR:("last_tag"/"tag_name"):"strand":"$4-last_end";"last_name"/"type_name";"last_dir"/"dir
									}	
								}
							
								last_end=$5
								last_name=type_name
								last_tag=tag_name
								last_dir=$7

		}' | sed 's/ /	/g' > $IGR_ACC_FILE
							

		cat $GENE_ACC_FILE >> $GENE_FILE
		
		cat $GENE_ACC_FILE $IGR_ACC_FILE | tr '\r' '\n' | awk -F"\t" -v REP=$REP '{$3="feature";$1=REP; print $0}' | sed 's/ /	/g' | sed 's/full_tag=transcript/full_tag=ncRNA/g' | sed 's/full_tag=tmRNA/full_tag=ncRNA/g' | sed 's/type=tmRNA/type=ncRNA/g' | sort -k4n >> $ALL_FILE
		
		REP_PARSED=`echo $REP | cut -d"." -f1`
		typeset -i NUM_FNA_REP
		NUM_FNA_REP=`cat $TEMP_FNA | grep ">" | grep -w $REP_PARSED | wc -l`

		if [ $NUM_FNA_REP == 1 ];then
			if [ $TEMP_FNA == $PARSED_FNA ];then
				echo "$TEMP_FNA and $PARSED_FNA are the same file"
			fi
			echo "fna_extract $TEMP_FNA $REP_PARSED $REP >> $PARSED_FNA"
			fna_extract $TEMP_FNA $REP_PARSED $REP >> $PARSED_FNA
		fi
	
		if [ $NUM_FNA_REP == 0 ];then
	
			echo "No fasta header with $REP_PARSED found in $TEMP_FNA!"
			cat $TEMP_FNA | grep ">" 
			exit
		fi
		if [ $NUM_FNA_REP -gt 1 ];then
			echo "$REP_PARSED has more than one fasta header!"
			cat $TEMP_FNA | grep ">" 
			exit
		fi

done

echo "Below are a list of feature types in "$ALL_FILE
while read line
do
	gff_info_parse "$line" 9 "type"
done < $ALL_FILE | sort | uniq -c | sort 
