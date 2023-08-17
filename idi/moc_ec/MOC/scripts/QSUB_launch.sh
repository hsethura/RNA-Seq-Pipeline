#!/bin/sh

echo  source /broad/IDP-Dx_storage/MOC/scripts/bash_header > qsub_launch.txt
echo  qsub -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat6 $@ >> qsub_launch.txt