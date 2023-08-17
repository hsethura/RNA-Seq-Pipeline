#/bin/sh

DIR=$1


ALL_FILES=`ls $DIR/* | grep count.txt | sort `


for FILE in $ALL_FILES
do	
	ROOT=`basename $FILE _read_count.txt`
	export TMPDIR=/broad/hptmp/livny/
	Total_cDNA=`cat $FILE | grep cdna | awk -v total=0 '{
			
															total=total+$2
															print total
															}' | tail -1`
	Total_hrRNA=`cat $FILE | grep rRNA | grep -v Mt_rRNA | awk -v total=0 '{
			
															total=total+$2
															print total
															}' | tail -1`
	Total_mtRNA=`cat $FILE | grep Mt_rRNA | awk -v total=0 '{
			
															total=total+$2
															print total
															}' | tail -1`

	pcnt_rRNA=`echo $Total_hrRNA $Total_mtRNA $Total_cDNA | awk '{print ($1+$2)/$3*100}'`
	pcnt_hrRNA=`echo $Total_hrRNA $Total_cDNA | awk '{print $1/$2*100}'`
	pcnt_mtRNA=`echo $Total_mtRNA $Total_cDNA | awk '{print $1/$2*100}'`
	
	echo $ROOT $pcnt_rRNA $pcnt_hrRNA $pcnt_mtRNA
 
done