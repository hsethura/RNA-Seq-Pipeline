# UGE_SUBMISSIONS

## UGER_cmd_batch_processor.py

UGER_cmd_batch_processor.py simplifies running large numbers of commands on UGE by batching commands into groups and tracking success/failure for each command within a group.  If the job is interrupted or failures occur, the UGER_cmd_batch_processor.py command can be rerun and only jobs that either failed or did not complete in the earlier execution will be resubmitted.

### Running UGER_cmd_batch_processor.py

UGER_cmd_batch_processor.py requires as input a file containing all the commands to be executed. Each command should be written as a unix-style shell command, and can direct stdout and stderr to any files it would like.  By default, the stdout and stderr for each command will be captured seprately within a job tracking directory created by the script itself.

For example, a 'commands.txt' file might include the following

    echo hello > hello.txt
    echo world > world.txt
    ls > files.list

Of course, your commands might do something more useful, like search a genome for mutations.  ;-)

>Note, it is best to have each command be responsible for writing its output to a designated location (ie. its own output file).

Execute these commands on UGE like so:

    UGER_cmd_batch_processor.py --cmds_file commands.txt \
                                --batch_size 2 \
                                --queue short \
                                --memory 1 \
                                --job_name my_job_name \
                                --tracking_dir tmp.tracking \
                                --project_name broad


The parameters are described as follows:

* --cmds_file: your file containing the commands to execute
* --batch_size:  the number of commands to bundle together and execute serially on the same grid node.
* --queue: the UGE queue to submit to.
* --memory: number of GB of RAM to reserve.
* --job_name: the name to give your job, as will be seen in a qstat report.
* --tracking_dir: the path to a directory that UGER_cmd_batch_processor.py will create and use for job management.
* --project_name: put your UGE project name (required to set your priority level).  If you don't have one, use 'broad'.

Additional options available include:

* --bash_header_text: this is a good place to put any dotkit 'reuse' statements or set up any special environament. ex. --bash_header_text "source /broad/software/scripts/useuse; reuse Java-1.8; reuse GCC-4.9;"  Note, if this is extensive, you can use the --bash_header parameter and point to a file that contains the full bash script header to include.
* --bash_header: point to a file containing the bash header to be executed prior to running your command, setting up your environment (see above alternative: --bash_header_text)
* --num_cores: number of threads required for processes that use multithreading
* --max_concurrent: maximum number of task array jobs to run concurrently.


### Job monitoring:

Use 'qstat' to monitor your jobs.  From the qstat report, it's straightforward to see how many of your jobs are running and how many remain queued.

     job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID 
     ------------------------------------------------------------------------------------------------------------------------------------------------
     4024350 0.34872 batch_runn bhaas        r     11/23/2016 11:10:51 long@uger-d001.broadinstitute.                                    1 1
     4024350 0.34872 batch_runn bhaas        r     11/23/2016 11:11:10 long@ugerbm-c001.broadinstitut                                    1 2
     4024350 0.32184 batch_runn bhaas        qw    11/23/2016 11:10:43                                                                   1 3-6084:1

Here, 2 of my jobs are running, and the remaining 3 through 6084 are pending.  Soon, dozens will be running in parallel.


### Tracking directory explained

The --tracking_dir doesn't normally require much attention, but it is useful to know how it is being used, particularly when needing to troubleshoot failed jobs.

The structure of the tracking directory is:

     -rwxrwx--- 1 bhaas broad 715 Nov 23 11:16 batch_runner.sh*
     -rw-rw-r-- 1 bhaas broad   0 Nov 23 11:16 cmds.FAILURES
     -rw-rw-r-- 1 bhaas broad   0 Nov 23 11:16 cmds.UNKNOWN
     drwxrwxr-x 2 bhaas broad 231 Nov 23 11:16 cmds/
     drwxrwxr-x 2 bhaas broad   0 Nov 23 11:16 ret/
     -rw-rw-r-- 1 bhaas broad  35 Nov 23 11:16 cmds.cached_success


The cmds/ directory contains a shell script for each bundle of commands that gets executed on a single node.  Any stderr or stdout (not redirected in a command-specific way according to your commands.txt file) will be written to the corresponding .err and .out file in this directory.  The number of cmd files will correspond to the total number of 'batches' generated and is determined by the total number of commands in your commands.txt file and the --batch_size setting.

The ret/ directory contains the exit value for each of the individual commands being executed. Note, the content of this directory is dynamic, and the UGER_cmd_batch_processor.py script routinely polls these file to collect them as they appear.

The commmands that complete successfully (as determined by an exit 0 code from running that single command) are written to the 'cmds.cached_success' file. Any failed commands (exiting non-zero) are written to the 'cmds.FAILURES' file.  If, when the initial UGE qsub array job command completes, there are commands that are not accounted for (missing .ret files), these are written to the 'cmds.UNKNOWN' file; these 'missing' jobs can occur if a batch of commands times-out on the node or there's some yet to be explained system glitch.

The 'batch_runner.sh' script will contain any settings provided by the --batch_header or --batch_header_text parameters, and this script is executed by UGE to execute the batched cmds in the cmds/ directory, based on the SGE_TASK_ID set by the job array.



### Resuming failed or unknown jobs

If anything should go awry and your UGER_cmd_batch_processor.py fails, you can rerun your original command, and it should only attempt to rerun the commands found in the cmds.FAILURES or cmds.UNKNOWN files.  Earlier outputs in the ret/ and cmds/ folder will be archived in the tracking_directory for later study.


