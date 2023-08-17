#!/bin/sh

i=0

while [ $i -le 100 ]
do
	SAM="/broad/hptmp/MOC/Acb_R_Gent_-_030m_rep1_NC_017847_pe.sam"
	COORDS="/broad/hptmp/Dx/test//Acb_R_Gent_-_030m_rep1_NC_017847_pe_coords_"$i".txt"
	METS="/broad/hptmp/Dx/test/Acb_R_Gent_-_030m_rep1_NC_017847_pe.mets_"$i".txt"
	ACC="/broad/hptmp/Dx///Acb_R_Gent_-_030m_rep1_NC_017847_pe_ACC.txt"
	QSUB_FILE="/broad/hptmp/Dx/test/qsub_script_"$i".txt"

	echo $COORDS

	echo "source /broad/IDP-Dx_work/nirmalya/bash_header" > $QSUB_FILE

	echo "/broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts//SAM_PARSE_old $SAM $COORDS $METS $ACC N " >>  $QSUB_FILE
	
	qsub -l h_rt=3:00:00 $QSUB_FILE

	i=`expr $i + 1`
done

exit

i=0

while [ $i -le 100 ]
do
	p=`expr $i + 1`
	while [ $p -le 100 ]
	do
	
		COORDS1="/broad/hptmp/Dx/test//Acb_R_Gent_-_030m_rep1_NC_017847_pe_coords_"$i".txt"
		COORDS2="/broad/hptmp/Dx/test//Acb_R_Gent_-_030m_rep1_NC_017847_pe_coords_"$p".txt"
		
		echo "diff $COORDS1 $COORDS2"
		diff $COORDS1 $COORDS2

		p=`expr $p + 1`
	done
	i=`expr $i + 1`
done