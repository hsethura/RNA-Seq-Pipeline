#!/bin/sh

NCBI_PATH=$1"/"
shift
OUT_DIR=$1"/"
shift
TEMP_DIR=$1"/"
shift
ADD5=$1
shift
ADD3=$1
shift

############################

TR_GFF_PARSE ()
{
	
	INFILE=$1
	GENE_KEY=$2
	PRODUCT_KEY=$3
	ACC=$4
	OUTFILE=$5
	APPEND=$6
	
	if [ $APPEND == "N" ];then
		
		echo "" | sed 1d > $OUTFILE
	
	fi
	
	cat $INFILE | tr '\r' '\n' | awk -v GENE_KEY=$GENE_KEY -v ACC=$ACC '{
							
						if($1 == ACC) 
						{	
							if($3 == GENE_KEY) 
							{
								y=split($9, ar, ";")
								
								for(i=1; i < y+1; i++)
								{
									if(ar[i] ~ /'$PRODUCT_KEY'/)
									{
										t=split(ar[i], pr, "=")
										PRODUCT=pr[2]
									}
									if(ar[i] ~ /gbkey=/)
									{
										t=split(ar[i], pr, "=")
										KEY=pr[2]
									}
									if(ar[i] ~ /ID=/)
									{
										t=split(ar[i], pr, "=")
										ID=pr[2]
									}
								}
								if(KEY=="")
									print $0";full_tag="ID":"PRODUCT":"ID":"$7":"$5-$4
								else
									print $0";full_tag="KEY":"PRODUCT":"ID":"$7":"$5-$4

							}
							
						}
					}' | sed 's/\%//g' >> $OUTFILE

 }

###########################

fna_extract ()
{

	FNA_FILE=$1
	ACC_PARSED=$2
	GFF_ACC=$3
	
	typeset -i NUM_MATCH
	NUM_MATCH=`	awk -v flag="N" -v GFF_ACC=$GFF_ACC -v ACC_PARSED=$ACC_PARSED '{
			if($1 ~ /'$ACC_PARSED'/ && $1 ~ />/ )
			{	
				print ">"GFF_ACC
			}
		}' $FNA_FILE | wc -l`

	
	cat $FNA_FILE | tr '\r' '\n' | awk -v flag="N" -v GFF_ACC=$GFF_ACC -v ACC_PARSED=$ACC_PARSED -v NUM_MATCH=$NUM_MATCH '{
	
					
					if($1 ~ />/ && flag=="Y")
						flag="N"
					if(flag=="Y" && $1 !~ />/)
 						print $1
					if ($1 ~ />/ )
					{
					
						if(NUM_MATCH==1)
						{
					
							if($1 ~ /'$ACC_PARSED'/)
							{	
								print ">"GFF_ACC
								flag="Y"
							}
						}
						if(NUM_MATCH > 1)
						{
							if($1 == ">"ACC_PARSED)
							{	
								print ">"GFF_ACC
								flag="Y"
							}
						}
					}
	}' | sed '/^$/d' 

}

############################

OUT_TAG=`echo $@ | sed 's/ /:/g'`
COMBINED_ALL_GFF=$OUT_DIR$OUT_TAG"_ALL.gff"
COMBINED_GENE_GFF=$OUT_DIR$OUT_TAG"_GENES.gff"
COMBINED_FNA=$OUT_DIR$OUT_TAG".fna"
COMBINED_VERS=$OUT_DIR$OUT_TAG"_fna_gff_version.txt"

echo "" | sed 1d > $COMBINED_ALL_GFF
echo "" | sed 1d > $COMBINED_GENE_GFF
echo "" | sed 1d > $COMBINED_FNA 
echo "" | sed 1d > $COMBINED_VERS

echo $OUT_TAG


for p in $@
do
	
	GFF_FILE=$NCBI_PATH"/"$p".gff"
	PARSED_GFF_FILE=$TEMP_DIR"/"$p"_parsed.gff"
	FNA_FILE=$NCBI_PATH"/"$p".fna"

	echo $GFF_FILE
	cat $GFF_FILE | sed 's/ /_/g' > $PARSED_GFF_FILE
		
	ALL_GFF_ACC=` awk '{if($1 ~ /FASTA/) exit; print $1}' $PARSED_GFF_FILE | sort | uniq | grep -v "#"  `
	
	echo $ALL_GFF_ACC
			
	for ACC in $ALL_GFF_ACC
	do	
		
		echo $ACC
		mkdir $TEMP_DIR"/"
								
		TEMP_CDS_FILE=$TEMP_DIR"/"$ACC"_acc_CDS_temp.gff"
		CDS_FILE=$TEMP_DIR"/"$ACC"_acc_CDS.gff"
		TRNA_FILE=$TEMP_DIR"/"$ACC"_acc_tRNA.gff"
		RRNA_FILE=$TEMP_DIR"/"$ACC"_acc_rRNA.gff"
		NCRNA_FILE=$TEMP_DIR"/"$ACC"_acc_ncRNA.gff"
		IGR_FILE=$TEMP_DIR"/"$ACC"_acc_IGR.gff"
		GENE_FILE=$OUT_DIR"/"$ACC"_acc_GENES.gff"
		ALL_FILE=$OUT_DIR"/"$ACC"_acc_ALL.gff"

		TR_GFF_PARSE $PARSED_GFF_FILE CDS product= $ACC $TEMP_CDS_FILE N
				
		echo $TEMP_CDS_FILE 
		echo $PARSED_GFF_FILE
		
		echo "Working on $TEMP_CDS_FILE for "$ACC" to generate "$CDS_FILE"..."
			
		while read line
		do  
			ACC_LINE=`echo $line | awk '{print $1}'`
			
			if [ $ACC == $ACC_LINE ];then	
			
				ID=`echo $line | awk '{	
											print "NONE"
											y=split($9, ar, ";")
											for(i=1; i < y+1; i++)
											{
												if(ar[i] ~ /Parent=/)
												{	
													split(ar[i], pr, "=")
													print pr[2]
													exit
												}
											}
										}' | tail -1`
			
				if [ $ID != 'NONE' ];then
			
					tag=`cat $PARSED_GFF_FILE | grep -w $ID | awk '{
															if($3 == "gene") 
															{	
																y=split($9, ar, ";")
																for(i=1; i < y+1; i++)
																{
																	if(ar[i] ~ /locus_tag=/)
																	{	
																		split(ar[i], pr, "=")
																		print pr[2]
																	}
																}
														
															}	
													}' | head -1`
			
			
					gene_name=`cat $PARSED_GFF_FILE | grep -w $ID | awk '{
															if(NR==1)
																print "-"
															if($3 == "gene") 
															{	
																y=split($9, ar, ";")
																for(i=1; i < y+1; i++)
																{
																	if(ar[i] ~ /gene=/)
																	{	
																		split(ar[i], pr, "=")
																		print pr[2]
																	}
																}
														
															}	
													}' | tail -1`
				else
			
					tag="-"
					gene_name="-"								
				fi
				echo $line | awk -v ADD5=$ADD5 -v ADD3=$ADD3 '{
											
										if($7=="+")
										{
											$4=$4-ADD5;$5=$5+ADD3
										}
										if($7=="-")
										{
											$4=$4-ADD3;$5=$5+ADD5
										}
										if($4 < 0)
											$4=0	
										print $0
									}' |  sed 's/full_tag=CDS/full_tag=CDS:'$tag':'$gene_name'/g' 
			
			fi
		done < $TEMP_CDS_FILE 	> $CDS_FILE 
				
		echo "Running TR_GFF_PARSE for "$ACC"..."

		TR_GFF_PARSE $PARSED_GFF_FILE tRNA product= $ACC $TRNA_FILE N
		TR_GFF_PARSE $PARSED_GFF_FILE rRNA product= $ACC $RRNA_FILE N
		TR_GFF_PARSE $PARSED_GFF_FILE ncRNA gene= $ACC $NCRNA_FILE N
		TR_GFF_PARSE $PARSED_GFF_FILE tmRNA gene= $ACC $NCRNA_FILE Y
		TR_GFF_PARSE $PARSED_GFF_FILE transcript product= $ACC $NCRNA_FILE Y
		
		echo $CDS_FILE
		echo $TRNA_FILE
		echo $RRNA_FILE
		echo $NCRNA_FILE
		echo $ALL_FILE
		echo $GENE_FILE		
		
		cat $TRNA_FILE $RRNA_FILE $CDS_FILE $NCRNA_FILE | sort -k4n | awk  '{type=$3;split($1, ar, ".");$1=ar[1];print $0";type="type}' | sed 's/ /	/g' > $GENE_FILE
		
		
		echo "" | sed 1d > $IGR_FILE
		
		cat $GENE_FILE | awk -v ADD5=$ADD5 -v ADD3=$ADD3 -v last_end=0 -v last_name=="ORIGIN" -v last_tag="ORIGIN" -v last_dir="|" '{
		
									$4=$4-1;$5=$5+1
									if($4 < 0)
										$4=0								
									
									
									y=split($9, ar, ";")
									type_field=ar[y]
									tag_field=ar[y-1]
									split(tag_field, qr, "=")
									tag_name=qr[2]							
									split(type_field, pr, "=")
									type_name=pr[2]							
									dir=$7
									strand="+"
									
									if(NR >= 1)
									{
									 
										$3="IGR"
										if(last_end < $4)
										{	
											gsub(":", ",", last_tag)
											gsub(":", ",", tag_name)

											print $1, $2, $3, last_end, $4, $6, strand, $8, last_tag"/"tag_name"):"strand":"$4-last_end";"last_name"/"type_name"_"last_dir"/"dir
										}	
									}
									
									last_end=$5
									last_name=type_field
									last_tag=tag_field
									last_dir=$7
		
								}' | sed 's/ /	/g' | sed 's/full_tag=/full_tag=IGR:(/g' |  sed 's/type=/type=IGR:/g' > $IGR_FILE
								
		echo $IGR_FILE
		
		cat $PARSED_GFF_FILE | grep "#" > $ALL_FILE
		
			
		cat $GENE_FILE $IGR_FILE | tr '\r' '\n' | awk -v ACC=$ACC '{$3="feature";$1=ACC; print $0}' | sed 's/ /	/g' | sed 's/full_tag=transcript/full_tag=ncRNA/g' | sed 's/full_tag=tmRNA/full_tag=ncRNA/g' | sed 's/type=tmRNA/type=ncRNA/g' | sort -k4n >> $ALL_FILE
		
		cat $ALL_FILE | tr '\r' '\n' | sed '/^$/d' >> $COMBINED_ALL_GFF
		cat $GENE_FILE >> $COMBINED_GENE_GFF
		
		ACC_PARSED=`echo $ACC | cut -d"." -f1`
		typeset -i NUM_FNA_ACC
		NUM_FNA_ACC=`cat $FNA_FILE | grep ">" | grep -w $ACC_PARSED | wc -l`
		
		if [ $NUM_FNA_ACC == 1 ];then
		
			echo "fna_extract $FNA_FILE $ACC_PARSED $ACC >> $COMBINED_FNA"
			fna_extract $FNA_FILE $ACC_PARSED $ACC >> $COMBINED_FNA
		fi
		
		if [ $NUM_FNA_ACC == 0 ];then
		
			echo "No fasta header with $ACC_PARSED found!"
			cat $FNA_FILE | grep ">" 
			exit
		fi
		if [ $NUM_FNA_ACC -gt 1 ];then
			echo "$ACC_PARSED has more than one fasta header!"
			cat $FNA_FILE | grep ">" 
			exit
		fi

		# making record of acc version in $COMBINED_VERS	
		cat $PARSED_GFF_FILE | grep "##sequence-region" >> $COMBINED_VERS	
		cat $FNA_FILE | grep ">" >> $COMBINED_VERS	
	
		echo $COMBINED_ALL_GFF
	done

done

cp $COMBINED_ALL_GFF $TEMP_DIR
cp $COMBINED_FNA $TEMP_DIR
cp $COMBINED_VERS $OUT_DIR

echo $COMBINED_ALL_GFF
echo $COMBINED_FNA
echo $COMBINED_VERS

exit
