#!/bin/sh


FILE=$1


gene_find ()
{

	FILE=$1
	TAG=$2
	TYPE=$3
	OUTFILE=$4
	
	cat $FILE | head -1 > $OUTFILE
	cat $FILE | grep $TAG | grep "^"$TYPE":" | awk '{
						split($1, ar, ":")
						split(ar[4], qr, "_")
 	#					$1=qr[1]
						print $0
						
				}' >> $OUTFILE
}

total ()
{

	FILE=$1
	ID=$2

	cat $FILE | awk -v ID=$ID -v t=1 -v total=0 '{										
										
										if(NR==1)
										{
											for(i=2; i < NF+1; i++)
												if($i==ID)
													t=i
										}
										if(NR!=1)
											total=total+$t
											
										print total
									
								}' | tail -1
}




gene_find $FILE 5S rRNA 5S.txt 
gene_find $FILE 16S rRNA 16S.txt
gene_find $FILE 23S rRNA 23S.txt 
gene_find $FILE _ CDS CDS.txt 
gene_find $FILE ssrA ncRNA tmRNA.txt 
gene_find $FILE _ ncRNA ncRNA.txt 



ALL_IDS=`cat 16S.txt | head -1 | awk '{
										for(i=2; i < NF+1; i++)
											print $i
									}'`


#echo $ALL_IDS

for ID in $ALL_IDS
do
	
	printf "%s\t" $ID
	total_all=`total $FILE $ID`
	total_5S=`total 5S.txt $ID`
	total_16S=`total 16S.txt $ID`
	total_23S=`total 23S.txt $ID`
	total_tmRNA=`total tmRNA.txt $ID`
	total_ncRNA=`total ncRNA.txt $ID`
	total_CDS=`total CDS.txt $ID`
	ratio_5S_CDS=`echo $total_5S $total_CDS | awk '{print $1/$2}'`
	ratio_16S_CDS=`echo $total_16S $total_CDS | awk '{print $1/$2}'`
	ratio_23S_CDS=`echo $total_23S $total_CDS | awk '{print $1/$2}'`
	ratio_23S_16S=`echo $total_23S $total_16S | awk '{print $1/$2}'`
	ratio_tmRNA_CDS=`echo $total_tmRNA $total_CDS | awk '{print $1/$2}'`
	ratio_5S_total=`echo $total_5S $total_all | awk '{print $1/$2}'`
	ratio_16S_total=`echo $total_16S $total_all | awk '{print $1/$2}'`
	ratio_23S_total=`echo $total_23S $total_all | awk '{print $1/$2}'`
	ratio_tmRNA_total=`echo $total_tmRNA $total_all | awk '{print $1/$2}'`

	#echo $ratio_5S_total $ratio_16S_total $ratio_23S_total $ratio_tmRNA_total 
	echo $ratio_5S_CDS $ratio_16S_CDS $ratio_23S_CDS $ratio_tmRNA_CDS 
	#ratio_23S_16S
# 	printf "%s\t" $ID
#	echo $total_CDS $total_5S $total_16S $total_23S $total_tmRNA $total_ncRNA $total_all
done




# 
# 	
# cat $FILE | awk '{
# 					if(NR == 1 || ($1 ~ /^rRNA:/ && $1 ~ /23S_/))
# 					{	
# 						split($1, ar, ":")
# 						split(ar[4], qr, "_")
# 						$1=qr[1]
# 						print $0
# 					}
# 						
# 				}' > 23S_counts.txt			
# 
# 
# i=1
# printf "Gene\t"
# 
# while [ $i -lt $NF ];
# do
# 
# 	name=`cat 16S_counts.txt | head -1 | awk -v total=0 -v i=$i '{print $i}'`
# 
# 	printf "%s\t" $name
# 	
# 	i=`expr $i + 1`
# 
# done
# echo ""
# 
# 
# printf "16S\t"
# 
# i=1
# 
# while [ $i -lt $NF ];
# do
# 
# 	total=`cat 16S_counts.txt | sed 1d | awk -v total=0 -v i=$i '{total=total+$i; print total}' | tail -1`
# 	
# 	printf "%s\t" $total
# 	
# 	i=`expr $i + 1`
# 
# done
# 
# echo ""
# printf "23S\t"
# 
# i=1
# 
# while [ $i -lt $NF ];
# do
# 
# 	total=`cat 23S_counts.txt | sed 1d | awk -v total=0 -v i=$i '{total=total+$i; print total}' | tail -1`
# 	
# 	printf "%s\t" $total
# 	
# 	i=`expr $i + 1`
# 
# done
# 
# echo ""
