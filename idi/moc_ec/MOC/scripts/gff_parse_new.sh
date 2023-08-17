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
TAG_FILE=$1
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

ls -lrt $GFF_FILE
ls -lrt $FNA_FILE

# output files
ALL_FILE=$OUT_DIR"/"$ACC"_ALL.gff"
GENE_FILE=$OUT_DIR"/"$ACC"_GENES.gff"
PARSED_FNA=$OUT_DIR$ACC".fna"
TEMP_FNA=$TEMP_DIR"/"$ACC"_temp.fna"


cat $GFF_FILE | tr '\r' '\n' | grep -i "^#" > $ALL_FILE
 
cat $FNA_FILE | sed -e 's/(//g' -e 's/)//g' -e 's/|/ /g' > $TEMP_FNA

echo "" | sed 1d > $GENE_FILE
echo "" | sed 1d > $PARSED_FNA 


ALL_GFF_REPS=`cat $GFF_FILE | tr '\r' '\n' | sed 1d | sed -e 's/ /_/g' -e 's/(//g' -e 's/)//g' -e 's/|/ /g' | awk '{if($1 ~ /FASTA/) exit; print $1}' |  sed '/^$/d' | sort -r | uniq | grep -v "#" `

for REP in $ALL_GFF_REPS
do	
	echo "REP: " $REP
	
	TEMP_ACC_FILE=$TEMP_DIR"/"$REP"_rep_temp.gff"
	GFF_ACC_FILE=$TEMP_DIR"/"$REP"_rep.gff"
	GENE_ACC_FILE=$TEMP_DIR"/"$REP"_rep_GENES.gff"
	IGR_ACC_FILE=$TEMP_DIR"/"$REP"_rep_IGR.gff"
	
	cat $GFF_FILE | tr '\r' '\n' | grep -v "#" | sed -e 's/ /_/g' -e '/^$/d' -e 's/)//g' -e 's/(//g' -e 's/|//g' -e 's/Coding_gene/CDS/g' | sort -k4n  | awk -v REP=$REP '{if($3 !~ /region/ && $3 !~ /ource/ && $1 == REP) print $0}' > $TEMP_ACC_FILE


	echo $GFF_ACC_FILE

	./GFF_parse $TEMP_ACC_FILE > $GFF_ACC_FILE

	echo $GFF_ACC_FILE


	TAGS=`cat $TAG_FILE | awk '{printf "%s;", $1}'`

	echo $TAGS

	cat $GFF_ACC_FILE | awk -v TAGS=$TAGS '{

							for(i=1; i < NF;i++)
								printf "%s\t", $i
														
							num_nf=split($NF, ar, ";")
							num_tags=split(TAGS, pr, ";")
						
							printf "type=%s;", $3
							printf "len=%s;", $5-$4

							for(p=1; p < num_tags; p++)
							{
								value="-"
								for(i=1; i < num_nf+1; i++)
								{	
									split(ar[i], sr, "=")
										if(sr[1]==pr[p])
								
									if(value=="-")
										value=sr[2]
									else
									{
										if(value!=sr[2])
											value=value","sr[2]
									}
								}
								
								if(value != "-")
									printf "%s=%s:", pr[p], value

							}
							print ""
						}' > $GENE_ACC_FILE 


	echo $GENE_ACC_FILE
	
	
	echo "Parsing $GENE_ACC_FILE to generate $IGR_ACC_FILE..."

	echo "" | sed 1d > $IGR_ACC_FILE

	cat $GENE_ACC_FILE | awk -F"\t" -v last_end=1 -v last_name=="ORIGIN" -v last_tag="ORIGIN" -v last_dir="|" '{

								$4=$4-1;$5=$5+1
								if($4 < 0)
									$4=0								
							
							
								y=split($9, ar, ";")
								split(ar[1], pr, "=")
								type_name=pr[2]
								tag_name=ar[3]
								len=$4-last_end
								igr_type=last_name"/"type_name":"last_dir"/"dir
								igr_ann="("last_tag"/"tag_name")"
								
								dir=$7
								strand="+"
							
								if(NR >= 1)
								{
							 
									$3="IGR"
									if(last_end < $4)
									{	
										gsub(":", ",", last_tag)
										gsub(":", ",", tag_name)

										print $1, $2, $3, last_end, $4, $6, strand, $8, "type=IGR;len="len";IGR_type="igr_type";IGR_ann="igr_ann
									}	
								}
							
								last_end=$5
								last_name=type_name
								last_tag=tag_name
								last_dir=$7

		}' | sed 's/ /	/g' > $IGR_ACC_FILE
							

		
		echo $IGR_ACC_FILE
		
		cat $GENE_ACC_FILE >> $GENE_FILE
			
		cat $GENE_ACC_FILE $IGR_ACC_FILE | tr '\r' '\n' | awk -F"\t" -v REP=$REP '{$3="feature";$1=REP; print $0}' | sed 's/ /	/g' | sed 's/type=transcript/type=ncRNA/g' | sed 's/type=tmRNA/full_tag=ncRNA/g'| sort -k5n >> $ALL_FILE
						
		
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
