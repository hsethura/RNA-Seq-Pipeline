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

cd `dirname $0`
export PIDFILE=monitor_jobs.pid
set -m
nohup ./uge_monitor_loop.sh >> monitor_loop.log 2>&1 &
set +m

# Wait a couple of seconds for python to launch
sleep 2
if [ ! -f $PIDFILE ]; then
	echo ERROR: python pid file not found: $PIDFILE.
	exit 1
fi
python_pid=`cat $PIDFILE`
echo "monitor_jobs launched with PID $python_pid"

# Wait a little while for it to maybe fail before checking for process existence.
sleep 10

if ps -p $python_pid > /dev/null
then echo monitor_jobs appears to be running.
else echo ERROR: monitor_jobs process not found. >&2
     exit 1
fi
