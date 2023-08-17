#!/bin/sh

GFF_FILE=$1
TAG_FILE=$2

TEMP_GFF=`echo $GFF_FILE | sed 's/.gff/_temp.gff/g'`
COMB_GFF=`echo $GFF_FILE | sed 's/.gff/_comb.gff/g'`
PARSED_GFF=`echo $GFF_FILE | sed 's/.gff/_parsed.gff/g'`

cat $GFF_FILE | sed 's/ /_/g'  | grep -v "#" > $TEMP_GFF

echo $TEMP_GFF

./GFF_parse $TEMP_GFF > $PARSED_GFF

echo $PARSED_GFF

#$TAG_FILE

TAGS="locus_tag;old_locus_tag;product;gene"
cat $PARSED_GFF | awk -v TAGS=$TAGS '{

							for(i=1; i < NF;i++)
								printf "%s\t", $i
							
							
							
							num_nf=split($NF, ar, ";")
							num_tags=split(TAGS, pr, ";")
						
							#print $NF
							for(p=1; p < num_tags+1; p++)
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
								printf "%s=%s;", pr[p], value

							}
							print ""
						}' 