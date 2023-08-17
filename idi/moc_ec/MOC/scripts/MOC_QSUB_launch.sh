#!/bin/sh

echo source /broad/IDP-Dx_work/nirmalya/bash_header > qsub.txt
echo $@ >> qsub.txt

qsub -l h_rt=24:00:00 qsub.txt


