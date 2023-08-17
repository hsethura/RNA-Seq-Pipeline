#!/bin/sh


FILE=$1


rRNA_find ()
{

	FILE=$1
	tag=$2
	OUTFILE=$3
	
	cat $FILE | head -1 > $OUTFILE
	cat $FILE | grep $tag | awk '{
					if(NR == 1)
						print "\t"$0

					if(NR != 1 && ($1 ~ /^rRNA:/))
					{	
						split($1, ar, ":")
						split(ar[4], qr, "_")
						$1=qr[1]
						print $0
					}
						
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


rRNA_find $FILE 16S 16S.txt 
rRNA_find $FILE 23S 23S.txt 

ALL_IDS=`cat 16S.txt | head -1 | awk '{
										for(i=2; i < NF+1; i++)
											print $i
									}'`


#echo $ALL_IDS

for ID in $ALL_IDS
do
	
	printf "%s\t" $ID
	total_16S=`total 16S.txt $ID`
	total_23S=`total 23S.txt $ID`
	ratio=`echo $total_23S $total_16S | awk '{print $1/$2}'`
	echo $ratio
	#echo $total_16S"\t"$total_23S"\t"$ratio
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
