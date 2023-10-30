#!/bin/bash

../UGER_cmd_batch_processor.py --cmds_file my_commands.txt \
                               --batch_size 2 \
                               --queue short \
                               --memory 1 \
                               --job_name my_uger_test \
                               --bash_header_text "source /broad/software/scripts/useuse; reuse Java-1.8; reuse GCC-4.9;" \
                               --tracking_dir tmp.tracking \
                               --project_name broad

