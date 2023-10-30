#!/usr/bin/env python
# encoding: utf-8

# contributed by bhaas@broadinstitute.org

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging
import argparse
import time
import subprocess


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="Processes commands using task array in batches and tracks success/failure at individual job level",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--cmds_file", dest="cmds_file", type=str, default="", required=True,
                        help="file containing all commands to be executed")

    parser.add_argument("--batch_size", dest="batch_size", type=int, default=0, required=True,
                        help="Number of commands to put in a single batch")

    parser.add_argument("--tracking_dir", dest='tracking_dir', type=str, default="", required=True,
                        help="directory for job tracking purposes")

    parser.add_argument("--project_name", dest="project_name", type=str, default="", required=True,
                        help="project name (shows up in qstat)")
    
    parser.add_argument("--queue", dest="queue", type=str, default="", required=False,
                        help="UGER queue name")
    parser.add_argument("--job_name", dest="job_name", type=str, default="no_name", help="name for job as viewed in qstat")

    parser.add_argument("--memory", dest="memory", type=int, required=True,
                        help="memory required (in GB) ex. --memory 10")
    
    parser.add_argument("--bash_header", dest='bash_header', type=str, default="", required=False,
                        help="include any 'use' statements of env settings to be used at the top of the submission script")

    parser.add_argument("--bash_header_text", dest="bash_header_text", type=str, default="", required=False)

    parser.add_argument("--debug", required=False, action="store_true", default=False, help="debug mode")

    parser.add_argument("--reservations", required=False, action="store_true", default=False,
                        help="turn on UGER resource reservations")

    parser.add_argument("--max_concurrent", required=False, type=int, default=None,
                        help="maximum number of concurrent processes")

    parser.add_argument("--num_cores", required=False, type=int, default=1,
                        help="number of cores to request for multithreading")
    parser.add_argument('--run_time', dest = 'run_time', default = "24", help = 'Run time for each node in hours')

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)      


    cmds_file = args.cmds_file
    tracking_dir = args.tracking_dir
    batch_size = args.batch_size
    run_time = args.run_time
    
    tracking_dir = os.path.abspath(tracking_dir)

    # initialize job management
    cached_successes = {}
    cached_success_file = os.path.sep.join([tracking_dir, "cmds.cached_success"])
    cached_success_fh = None

    cmds_dir = os.path.sep.join([tracking_dir, 'cmds'])
    rets_dir = os.path.sep.join([tracking_dir, 'ret'])
    
    if not os.path.isdir(tracking_dir):
        os.makedirs(tracking_dir)
        cached_success_fh = open(cached_success_file, 'w', buffering=0)
        
    else:
        # resume from earlier processing
        cached_successes = parse_cached_successes_file(cached_success_file)
        cached_success_fh = open(cached_success_file, 'a', buffering=0)

        # deprecate the old cmds and rets dirs
        timestamp = time.strftime("%m-%d-%y-t%H:%M:%S")
        reloc_cmds_dir = os.path.sep.join([tracking_dir, "__" + timestamp + "_cmds"])
        reloc_rets_dir = os.path.sep.join([tracking_dir, "__" + timestamp + "_rets"])

        os.rename(cmds_dir, reloc_cmds_dir)
        os.rename(rets_dir, reloc_rets_dir)

    # build new workspace
    os.makedirs(cmds_dir)
    os.makedirs(rets_dir)

    # capture the list of commands, ignore those that are already cached as successes
    cmds_list = parse_commands_file(cmds_file, cached_successes)

    if len(cmds_list) == 0:
        logger.info("All commands are already completed")
        sys.exit(0)

        
    ## write the batch runner
    bash_header_text = args.bash_header_text
    if (args.bash_header):
        bash_header_text = open(args.bash_header).read()
    
    batch_runner_script = os.path.sep.join([tracking_dir, "batch_runner.sh"])
    write_batch_runner_script(batch_runner_script, cmds_dir, bash_header_text)

    # batch the commands
    batch_count = batch_commands(cmds_list, batch_size, cmds_dir, rets_dir)

    uger_cmd = construct_uger_command(batch_runner_script, batch_count, args, cmds_dir, run_time)

    logger.info("UGER CMD: {}".format(uger_cmd))
    
    cached_failures_file = os.path.sep.join([tracking_dir, "cmds.FAILURES"])
    cached_failures_fh = open(cached_failures_file, 'w', buffering=0)
    cached_unknown_fh = open(os.path.sep.join(([tracking_dir, "cmds.UNKNOWN"])), 'w', buffering=0)


    # run the jobs
    (success_count, failure_count, unknown_count) = execute_and_track_uger_cmd(uger_cmd,
                                                                               rets_dir,
                                                                               cmds_list,
                                                                               cached_success_fh,
                                                                               cached_failures_fh,
                                                                               cached_unknown_fh)
    
    cached_failures_fh.close()
    cached_success_fh.close()
    cached_unknown_fh.close()

    # report on success/failure
    logger.info("Successful commands executed: {}".format(success_count))

    if unknown_count:
        logger.info("Unknown status count: {}".format(unknown_count))

    if failure_count:
        logger.info("Failure count: {}".format(failure_count))
                

    if failure_count + unknown_count > 0:
        errmsg = "Error, not all commands succeeded."
        logger.critical(errmsg)
        raise RuntimeError(errmsg)
    else:
        logger.info("All commands succeeded")
        sys.exit(0)
        
    


def execute_and_track_uger_cmd(uger_cmd, rets_dir, cmds_list, cached_success_fh, cached_failures_fh, cached_unknown_fh,
                               poll_waiting_time_init = 15, poll_waiting_max_time = 60):


    """
    all the hard work is done here...  runs the task-array job and tracks success/failure for each of the sub commands executed on the grid.
    """


    exit_val_tracker = []
    for i in range(0,len(cmds_list)):
        exit_val_tracker.append(-1)
        
    s = subprocess.Popen(uger_cmd, shell=True)

    poll_waiting_time = poll_waiting_time_init
    ret = s.poll()
    while ret == None:
        logger.debug("waiting for: {}".format(poll_waiting_time))
        time.sleep(poll_waiting_time)
        logger.debug("done waiting\n")
        poll_waiting_time = min(poll_waiting_max_time, poll_waiting_time *2)
        
        uger_rets_dict = find_ret_vals(rets_dir)
        
        for cmd_idx, exit_val in uger_rets_dict.items():
            cmd_idx = int(cmd_idx)
            exit_val = int(exit_val)
            cmd = cmds_list[cmd_idx]
            exit_val_tracker[cmd_idx] = exit_val
            logger.debug("cmd_idx: {} exited with {}".format(cmd_idx,exit_val))
            if exit_val == 0:
                cached_success_fh.write(cmd + "\n")
            else:
                cached_failures_fh.write(cmd + "\n")

        ret = s.poll()

    ret = s.wait() # collect child

    if ret:
        logger.error("uger cmd: {} died with exit code: {}".format(uger_cmd, ret))
        # note, no need to raise exception here. continue to report counts of job statuses

    
    success_count = 0
    failure_count = 0
    unknown_count = 0

    for i in range(0,len(cmds_list)):

        exit_code = exit_val_tracker[i]
        
        if exit_code == 0:
            success_count += 1
        elif exit_code > 0:
            failure_count += 1
        else:
            unknown_count += 1
            cached_unknown_fh.write(cmds_list[i] + "\n")

    return(success_count, failure_count, unknown_count)


        

def find_ret_vals(dir):
    """
    captures the exit values from each of the sub-commands that were executed
    """


    ret_files = subprocess.check_output("find {} -type f -regex \".*ret\"".format(dir), shell=True)

    if not ret_files:
        return {}

    uger_rets_dict = {}
    
    ret_files_list = ret_files.split("\n")

    for ret_file in ret_files_list:
        logger.debug("found ret file: {}".format(ret_file))
        x = re.search("/(\d+)\.ret$", ret_file)
        if x:
            cmd_idx = x.group(1)
            ret_val = open(ret_file).read()
            ret_val = ret_val.strip()

            logger.debug("idx {} exited w/ ret_val: {}\n".format(cmd_idx, ret_val))
            
            uger_rets_dict[cmd_idx] = ret_val
            os.unlink(ret_file)

    return uger_rets_dict
    


def construct_uger_command(batch_runner_script, batch_count, args, cmds_dir, run_time):
    """
    constructs the main UGE command that runs in task-array mode
    """


    time_str = run_time + ":00:00"
    cmd = str ("qsub -V -cwd -b y -sync y " +
               #" -q " + args.queue +
               "  -l h_rt=" + time_str + 
               " -pe smp " + str(args.num_cores) + " -binding linear:" + str(args.num_cores) +
               #" -l os=RedHat6 " + 
               " -l h_vmem=" + str(args.memory) + "g " +
               " -P " + args.project_name +
               " -N " + args.job_name + 
               " -e " + os.path.dirname(cmds_dir) + ".qsub.err " +
               " -o " + os.path.dirname(cmds_dir) + ".qsub.out ")
    
    if args.reservations:
        cmd += " -R y "

    cmd += " -t 1-" + str(batch_count)

    if args.max_concurrent:
        cmd += " -tc " + str(args.max_concurrent)

    cmd += " " + batch_runner_script

    return cmd

    

def parse_cached_successes_file(cached_success_file):

    cached_successes = {}

    with open(cached_success_file, 'r') as f:
        for line in f:
            cmd = line.rstrip()
            cached_successes[cmd] = True

    return cached_successes
    

    

def parse_commands_file(cmds_file, cached_successes):

    cmds_list = []

    with open(cmds_file) as f:
        for line in f:
            cmd = line.rstrip()
            if re.search("\w", cmd) and not cached_successes.has_key(cmd):
                cmds_list.append(cmd)

    return cmds_list


def write_batch_runner_script(script_filename, cmds_dir, bash_header_text=None):

    ofh = open(script_filename, 'w')

    ofh.write("#!/bin/bash\n\n")
    if bash_header_text:
        ofh.write(bash_header_text + "\n\n")
    
    script_contents = """
if [ -z \"\$SGE_TASK_ID\" ]; then
    echo \"ERROR, env var SGE_TASK_ID is not set\"
    exit 1
fi

echo \"HOSTNAME ${HOSTNAME}, SGE_TASK_ID: ${SGE_TASK_ID}\" > __CMDs_DIR__/${SGE_TASK_ID}.cmd.err

cmd=\"__CMDs_DIR__/${SGE_TASK_ID}.cmd 1>__CMDs_DIR__/${SGE_TASK_ID}.cmd.out 2>>__CMDs_DIR__/${SGE_TASK_ID}.cmd.err\"

eval $cmd

exit 0
    """

    script_contents = script_contents.replace("__CMDs_DIR__", cmds_dir)


    
    ofh.write(script_contents)
    ofh.close()
    os.chmod(script_filename, 0770)

    

def batch_commands(cmds_list, batch_size, cmds_dir, rets_dir):

    batch_count = 0
    cmd_idx = -1

    
    while cmd_idx < len(cmds_list):

        batch_count += 1

        batch_cmd_file = os.path.sep.join([cmds_dir, str(batch_count) + ".cmd"])
        ofh = open(batch_cmd_file, 'w')
        ofh.write("#!/bin/bash\n\n")
        
        for i in range(0,batch_size):
            cmd_idx += 1
            if cmd_idx >= len(cmds_list):
                break
            cmd = cmds_list[cmd_idx]
            ofh.write("{}\n".format(cmd))
            ofh.write("echo $? > {}\n\n".format(os.path.sep.join([rets_dir,str(cmd_idx)+".ret"])))

        ofh.write("exit 0\n")
        ofh.close()
        os.chmod(batch_cmd_file, 0770)

    return batch_count


        

 
####################
 
if __name__ == "__main__":
    main()
