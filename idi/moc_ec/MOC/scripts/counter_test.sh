#!/bin/sh



i=1

while [ $i -lt 10 ]
do
	
	echo "source /broad/IDP-Dx_storage/MOC/scripts/bash_header" > counter_test.txt
	echo "/broad/IDP-Dx_work/nirmalya/pipeline/beta/PipelineMain.py -c /broad/IDP-Dx_work/nirmalya/pipeline/configs/config_default.yaml --key_path /broad/IDP-Dx_storage/MOC/Key_files/MOC595_aln_test_key.txt --raw_seq_path /broad/hptmp/vbetin/raw_seq_data --do_patho --min_resource --ADD3 30 --ADD5 20 --rm_rts_dup --MOC_id $i " >> counter_test.txt
	OUT_FILE="countertest_"$i"_o.txt"
	ERR_FILE="countertest_"$i"_e.txt"

	echo "qsub -o $OUT_FILE -e $ERR_FILE -l os=RedHat7 -l h_rt=24:00:00 "counter_test.txt""
	qsub -o $OUT_FILE -e $ERR_FILE -l os=RedHat7 -l h_rt=24:00:00 "counter_test.txt" 

	i=`expr $i + 1`


done