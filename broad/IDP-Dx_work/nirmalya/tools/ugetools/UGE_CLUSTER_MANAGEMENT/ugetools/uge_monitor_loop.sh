#!/usr/bin/env bash
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2017 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.
# 
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.

source /broad/software/scripts/useuse
unuse -q Python-2.7
use Python-2.7

cd `dirname $0`

exit_status=129
admin_email=alecw@broadinstitute.org

while (( $exit_status > 128 ))
do
    date
    # enable job control so that kill -INT can be sent to java process
    set -m
    set -x
    ./monitor_jobs.py --sleep-mins 130 --exclude-user root -v  -e uger-monitor@broadinstitute.org \
    --ignore-job-file ignore_job_file.monitor_jobs --exclude-user-file excluded_users.monitor_jobs  >> monitor_jobs.log 2>&1 &
    set +x

    python_pid=$!
    echo $python_pid > $PIDFILE

    set +m

    # This waits until the process terminates, and set $? to the exit status of that program
    wait $python_pid
    exit_status=$?

    if (( $exit_status > 128 ))
    then
     declare -i signal_number=$exit_status-128
     mail -s "Restarting monitor_jobs after terminated with signal" $admin_email <<EOF
Deployment Directory: `pwd`
Host: `hostname`
Terminated Process ID: $python_pid

monitor_jobs.py terminated with signal $signal_number
Restarting...
EOF
    fi
done

