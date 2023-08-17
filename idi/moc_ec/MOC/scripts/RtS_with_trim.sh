#!/bin/sh

MOC_ID=$1

KEY_FILE="/broad/IDP-Dx_storage/MOC/Key_files/"$MOC_ID"_key.txt"

# PROJ_ID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk '{print $1}' | sort | uniq`
# TEMP_DIR="/broad/hptmp/RNASeq_proj/MOC"

SPLIT_DIR="/broad/hptmp/RNASeq_proj/MOC/"$MOC_ID"/"*"/splitdir/"

ls -lrt $SPLIT_DIR

sh /broad/IDP-Dx_storage/MOC/scripts/RtS_Analysis_launch.sh $MOC_ID :--no_merge --no_align --no_count:

echo "sh /broad/IDP-Dx_storage/MOC/scripts/trim_fastq.sh 1 $SPLIT_DIR GATCTTTCA GTCGTCCGT GTGAATATT GCCAGGAGT"
sh /broad/IDP-Dx_storage/MOC/scripts/trim_fastq.sh 1 $SPLIT_DIR GATCTTTCA GTCGTCCGT GTGAATATT GCCAGGAGT

echo "sh /broad/IDP-Dx_storage/MOC/scripts/RtS_Analysis_launch.sh $MOC_ID :--no_split:"
sh /broad/IDP-Dx_storage/MOC/scripts/RtS_Analysis_launch.sh $MOC_ID :--no_split: