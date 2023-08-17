#!/bin/sh


BAM_INPUT=$1
REF_FILE=$2
BAMMET_OUTDIR=$3
PICARD_BIN=$4

mkdir -p $BAMMET_OUTDIR

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

ROOT=`echo $BAM_INPUT | rev | cut -d"/" -f1 | rev | sed 's/.pe.bam//g'`
OUT_MET=$BAMMET_OUTDIR"/"$ROOT"_AlignmentSummaryMetrics.txt"
OUT_QDIST=$BAMMET_OUTDIR"/"$ROOT"_MeanQualByCycle.txt"
OUT_QCHART=$BAMMET_OUTDIR"/"$ROOT"_MeanQualByCycle.pdf"
OUT_INSMET=$BAMMET_OUTDIR"/"$ROOT"_InsertSizeMets.txt"
OUT_INSHISTO=$BAMMET_OUTDIR"/"$ROOT"_InsertSizeHisto.pdf"

java -Xmx4g -jar $PICARD_BIN"/picard.jar" CollectAlignmentSummaryMetrics R=$REF_FILE I=$BAM_INPUT O=$OUT_MET VALIDATION_STRINGENCY=LENIENT
java -Xmx4g -jar $PICARD_BIN"/picard.jar" MeanQualityByCycle PF_READS_ONLY=false ALIGNED_READS_ONLY=false VALIDATION_STRINGENCY=LENIENT I=$BAM_INPUT O=$OUT_QDIST CHART=$OUT_QCHART
java -Xmx4g -jar $PICARD_BIN"/picard.jar" CollectInsertSizeMetrics I=$BAM_INPUT  O=$OUT_INSMET H=$OUT_INSHISTO M=0.5 VALIDATION_STRINGENCY=LENIENT




#sh $SCRIPTS_PATH"SATURATION.sh" $ALIGNED_SAM $ACC $ID $ACC_PATH

# 
# 
# java -jar $JAR_PATH CollectAlignmentSummaryMetrics R=$REF_FILE I=$ALIGNED_SAM O=/broad/hptmp/Dx/metrics.txt
# 
# java -jar $JAR_PATH MeanQualityByCycle I=$ALIGNED_SAM O=$BAMMET_OUTDIR"mean_qual_by_cycle.txt" CHART=$BAMMET_OUTDIR"mean_qual_by_cycle.pdf"
# 
# java -jar $JAR_PATH CollectQualityYieldMetrics I=$ALIGNED_SAM O=$BAMMET_OUTDIR"quality_yield_metrics.txt" 
# 
# java -jar $JAR_PATH CollectBaseDistributionByCycle CHART=$BAMMET_OUTDIR"collect_base_dist_by_cycle.pdf" I=$ALIGNED_SAM O=$BAMMET_OUTDIR"output.txt"
# 
# java -jar $JAR_PATH QualityScoreDistribution I=$ALIGNED_SAM O=$BAMMET_OUTDIR"qual_score_dist.txt" CHART=$BAMMET_OUTDIR"qual_score_dist.pdf"
# 
# java -jar $JAR_PATH SortSam I=$ALIGNED_SAM O=$ROOT"_coordSorted.sam" SORT_ORDER=coordinate
# 
# java -jar $JAR_PATH MarkDuplicates I=$ROOT"_coordSorted.sam" O=$ROOT"_dupMarked.sam" METRICS_FILE=$BAMMET_OUTDIR"dupMarked_metrics.txt"

