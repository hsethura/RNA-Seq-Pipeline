#!/bin/sh


PROJ_ID=$1
BAM_PATH=$2
OUT_PATH=$3
REF_FILE=$4

BAM_DIR=$BAM_PATH"/"$PROJ_ID
OUT_DIR=$OUT_PATH"/"$PROJ_ID"/BAM_metrics/"

mkdir -p $OUT_DIR

PICARD_BIN="/seq/software/picard/1.782/bin/"


mkdir $OUT_DIR
########  function for extracting options from command line
 
		extract_option ()
		{
			option_name=$1 
			shift
			default=$1
			shift 
			echo $@ | awk -v def=$default -v name=$option_name '{
		
				for(i=1; i < NF+1; i++)
				{
					if(i==1)			# intitalize to $default in case no option found
						print def
					if($i==name)
						print $(i+1)
				}
		
			}' | tail -1
		}

###############################  extract and checkpoint options  ##############################

options=$@

q=`extract_option -q short $options`


ALL_BAMS=`ls -lrt $BAM_DIR"/"*.bam | awk '{print $9}'`


for m in $ALL_BAMS
do

	ROOT=`echo $m | rev | cut -d"/" -f1 | rev | sed 's/.pe.bam//g'`
	OUT_MET=$OUT_DIR"/"$ROOT"_AlignmentSummaryMetrics.txt"
	OUT_QDIST=$OUT_DIR"/"$ROOT"_MeanQualByCycle.txt"
	OUT_QCHART=$OUT_DIR"/"$ROOT"_MeanQualByCycle.pdf"
	OUT_INSMET=$OUT_DIR"/"$ROOT"_InsertSizeMets.txt"
	OUT_INSHISTO=$OUT_DIR"/"$ROOT"_InsertSizeHisto.pdf"
	
	#java -jar $PICARD_BIN"/CollectAlignmentSummaryMetrics.jar" R=$REF_FILE I=$m O=$OUT_MET VALIDATION_STRINGENCY=LENIENT
	#java -jar $PICARD_BIN"/MeanQualityByCycle.jar" PF_READS_ONLY=false ALIGNED_READS_ONLY=false VALIDATION_STRINGENCY=LENIENT I=$m O=$OUT_QDIST CHART=$OUT_QCHART
	java -jar $PICARD_BIN"/CollectInsertSizeMetrics.jar" I=$m  O=$OUT_INSMET H=$OUT_INSHISTO M=0.5 VALIDATION_STRINGENCY=LENIENT 
	echo $OUT_INSMET
exit

done


exit


exit

#sh $SCRIPTS_PATH"SATURATION.sh" $ALIGNED_SAM $ACC $ID $ACC_PATH

exit


java -jar $JAR_PATH CollectAlignmentSummaryMetrics R=$REF_FILE I=$ALIGNED_SAM O=/broad/hptmp/Dx/metrics.txt

java -jar $JAR_PATH MeanQualityByCycle I=$ALIGNED_SAM O=$OUT_DIR"mean_qual_by_cycle.txt" CHART=$OUT_DIR"mean_qual_by_cycle.pdf"

java -jar $JAR_PATH CollectQualityYieldMetrics I=$ALIGNED_SAM O=$OUT_DIR"quality_yield_metrics.txt" 

java -jar $JAR_PATH CollectBaseDistributionByCycle CHART=$OUT_DIR"collect_base_dist_by_cycle.pdf" I=$ALIGNED_SAM O=$OUT_DIR"output.txt"

java -jar $JAR_PATH QualityScoreDistribution I=$ALIGNED_SAM O=$OUT_DIR"qual_score_dist.txt" CHART=$OUT_DIR"qual_score_dist.pdf"

java -jar $JAR_PATH SortSam I=$ALIGNED_SAM O=$ROOT"_coordSorted.sam" SORT_ORDER=coordinate

java -jar $JAR_PATH MarkDuplicates I=$ROOT"_coordSorted.sam" O=$ROOT"_dupMarked.sam" METRICS_FILE=$OUT_DIR"dupMarked_metrics.txt"

