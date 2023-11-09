#!/bin/bash

gzip -d -c /broad/hptmp/RNASeq_proj/MOC///hsethura//MOCP-0112//JJ5_17_23/mergedir//0PC_3_R2.fastq.gz > Y/adapter_analysis/0PC_3_R2.fastq
echo $? > /home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/ret/70.ret

awk 'NR%4==1{printf ">%s\n", substr($0,2)} NR%4==2{print}' Y/adapter_analysis/0PC_3_R2.fastq > Y/adapter_analysis/0PC_3_R2.fasta
echo $? > /home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/ret/71.ret

blastn -db /home/unix/hsethura/RNA-Seq-Pipeline/idi/moc_ec/MOC/RtS_pipeline/read_insertions_in_adapter/rts_adapters.fasta -query Y/adapter_analysis/0PC_3_R2.fasta -out Y/adapter_analysis/0PC_3_R2_blast_op_ws_10.csv -outfmt 10 -ungapped  -word_size 10 -strand plus -num_threads 10
echo $? > /home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/ret/72.ret

python3 /home/unix/hsethura/RNA-Seq-Pipeline/idi/moc_ec/MOC/RtS_pipeline/read_insertions_in_adapter/process_blast_output.py -csv Y/adapter_analysis/0PC_3_R2_blast_op_ws_10.csv -fasta Y/adapter_analysis/0PC_3_R2.fasta -lc_method RtS-TS
echo $? > /home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/ret/73.ret

rm Y/adapter_analysis/0PC_3_R2_blast_op_ws_10.csv
echo $? > /home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/ret/74.ret

rm Y/adapter_analysis/0PC_3_R2.fastq
echo $? > /home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/ret/75.ret

rm Y/adapter_analysis/0PC_3_R2.fasta
echo $? > /home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/ret/76.ret

exit 0
