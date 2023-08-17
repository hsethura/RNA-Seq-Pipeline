
#!/bin/sh

echo source /broad/IDP-Dx_work/nirmalya/bash_header > qsub.txt
echo /broad/IDP-Dx_work/nirmalya/pipeline/beta/PipelineMain.py $@ >> qsub.txt

qsub -l h_rt=4:00:00 qsub.txt



