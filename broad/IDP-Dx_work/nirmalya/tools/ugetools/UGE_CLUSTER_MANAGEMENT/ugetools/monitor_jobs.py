#!/usr/bin/env python
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2016 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.
"""
Daemon program that examines Grid Engine jobs in order to determine if any hosts are unhealthy.

A list of jobs for all users (by default) is captured, and the following conditions are reported:

- Host for a running job does not respond to ping
- Host for a running job has a queue in unreachable stage
- Job has been in transfer state too long
- Job has been in deletion stage too long
- Job is using more CPU than it has requested
- Job does appears to be cpu-starved; i.e. cpu-time/run-time is small
  Note that jobs can appear cpu-starved for legitimate reasons.  Therefore a host must have multiple cpu-starved jobs
  to be reported, and if a job array has cpu-starved jobs on multiple hosts, it is not considered to be a problem.

See the options below to control the thresholds for these conditions.

Note that files containing lists of users, hosts, jobs are re-read every time status is checked, so these files
can be changed without having to restart the daemon.

Note that files containing lists may have blank lines, and comments introduced with # and continuing to EOL.
"""
import argparse
import collections
import itertools
import logging
import signal
import subprocess
import sys
import time
import traceback
import email_notify

import uge_functions


def read_list_file(path):
    """
    Parse a file containing one element per line, optionally with line comments starting with '#'
    :param path: file to be read
    :return: list of elements, with comments and newlines removed
    """
    with open(path, "r") as fIn:
        ret = [line.rstrip("\n").split("#", 1)[0].strip() for line in fIn.readlines()]

    return [line for line in ret if len(line) > 0]


def combine_option_lists(list1, list2, default_list=None):
    if default_list is None:
        default_list = []
    if list1 is None:
        if list2 is None:
            if default_list is None:
                return []
            else:
                return default_list
        else:
            return list2
    elif list2 is None:
        return list1
    else:
        return list1 + list2


if sys.platform == "darwin":
    wait_flag = "W"
else:
    wait_flag = "w"

dev_null = open("/dev/null", "w")


def ping(host):
    """
    :arg host: host to be pinged
    :return: True if host responded to ping
    """
    # Send 1 packet and wait up to 3 seconds for a response.
    return subprocess.call(["ping", "-c", "1", "-" + wait_flag, "3", host],
                           stdout=dev_null, stderr=subprocess.STDOUT) == 0


def ping_hosts(hosts_to_ping, dct_no_ping_hosts):
    """
    :param hosts_to_ping: set of hosts to be pinged
    :param dct_no_ping_hosts: key: host; value: first time host failed ping.  Updated as appropriate
    :return: list of hosts that newly failed ping
    """
    ret = []
    for host in hosts_to_ping:
        if ping(host):
            if host in dct_no_ping_hosts:
                del dct_no_ping_hosts[host]
        elif host not in dct_no_ping_hosts:
            ret.append(host)
            dct_no_ping_hosts[host] = time.time()

    return ret


def check_host_queue_state(hosts, dct_bad_state_hosts):
    """
    Query queues on given hosts
    :param hosts: set of hosts to be queried
    :param dct_bad_state_hosts: key: host, value: first time host had bad state.  Updated as appropriate
    :return: list of hosts that newly have bad state
    """
    ret = []
    current_bad_state_hosts = uge_functions.get_hosts_with_unreachable_queues(hosts)
    for host in hosts:
        if host in current_bad_state_hosts:
            if host not in dct_bad_state_hosts:
                ret.append(host)
                dct_bad_state_hosts[host] = time.time()
        elif host in dct_bad_state_hosts:
            # problem has cleared
            del dct_bad_state_hosts[host]
    return ret


def is_transfer_state(state):
    return state.find("t") != -1


def is_delete_state(state):
    return state.find("d") != -1


def is_running_state(state):
    return state == "r"


ProblemState = collections.namedtuple("ProblemState", ["host", "state", "first_seen_time", "reported_time", "user",
                                                       "start_time"])


def update_strange_jobs(dct_strange_state_jobs, job, state, job_details=None):
    prev_problem_state = dct_strange_state_jobs.get(job.job_and_task)
    if prev_problem_state is None or prev_problem_state.state != state:
        # if state has changed, treat this as the old problem being cleared
        if job_details is not None:
            start_time_secs = job_details.start_time_secs
        else:
            start_time_secs = None
        dct_strange_state_jobs[job.job_and_task] = \
            ProblemState(job.host, state, time.time(), None, job.user, start_time_secs)


def cpu_starved(job_detail, wallclock_threshold_minutes, cpu_fraction):
    if job_detail.start_time_secs is None:
        return None
    wallclock_secs = time.time() - job_detail.start_time_secs
    if wallclock_secs < wallclock_threshold_minutes * 60:
        # job has not been running long enough
        return None
    elif job_detail.cpu_secs is None:
        logging.log_message("Could not get cpu time for ", job_detail.job_and_task)
        return 0
    else:
        fraction = job_detail.cpu_secs / float(wallclock_secs)
        if fraction <= cpu_fraction:
            return fraction
        else:
            return None


def cpu_hog(job_detail, hog_threshold):
    if job_detail.start_time_secs is None:
        return None
    wallclock_secs = time.time() - job_detail.start_time_secs
    if job_detail.cpu_secs is None:
        logging.log_message("Could not get cpu time for ", job_detail.job_and_task)
        return 0
    elif job_detail.slots is None:
        logging.log_message("Could not get slots for ", job_detail)
    else:
        hog_amount = (job_detail.cpu_secs / float(wallclock_secs)) - job_detail.slots
        if hog_amount >= hog_threshold:
            return hog_amount
        else:
            return None


def remove_strange_job(dct_strange_state_jobs, strange_job):
    if strange_job in dct_strange_state_jobs:
        del dct_strange_state_jobs[strange_job]


def count_starved_jobs_by_host(dct_strange_state_jobs,
                               job_per_host_verbosity_threshold=None,
                               host_per_job_verbosity_threshold=None):
    # Count the unique jobs on a host that are cpu-starved.
    # Multiple tasks for the same job on the same host count as 1.
    # The idea is that a job could look cpu-starved if it is waiting for something, e.g. downloading something
    # over the network.  If a host is in trouble, there should be multiple jobs that appear cpu-starved.
    # Also, count the number of hosts on which one of the tasks for a job is cpu-starved.  If there are several,
    # it is probably the job and not a host in trouble.
    # returns tuple(dict(host=>count of unique starved jobs), dict(job_id=>count of unique hosts on which a task for that job is starved)

    dct_host_starved_jobs = {}  # key: host; value: set of job IDs (not job_and_task)
    dct_starved_job_hosts = {}  # key: job ID, value: set of hosts on which that job has a starved task
    for job_and_task, problem_state in dct_strange_state_jobs.items():
        if problem_state.state.startswith("cpu-starved"):
            dct_host_starved_jobs.setdefault(problem_state.host, set()).add(job_and_task.job_id)
            dct_starved_job_hosts.setdefault(job_and_task.job_id, set()).add(problem_state.host)

    if job_per_host_verbosity_threshold is not None:
        for host, jobs in dct_host_starved_jobs.iteritems():
            if len(jobs) >= job_per_host_verbosity_threshold:
                logging.log_message(host + " has cpu-starved jobs " + ", ".join(jobs))

    if host_per_job_verbosity_threshold is not None:
        for job, hosts in dct_starved_job_hosts.iteritems():
            if len(hosts) >= host_per_job_verbosity_threshold:
                logging.log_message(job + " is cpu-starved on hosts " + ", ".join(hosts))

    return (dict([(host, len(jobs)) for host, jobs in dct_host_starved_jobs.items()]),
            dict([(job_id, len(hosts)) for job_id, hosts in dct_starved_job_hosts.items()]))


def format_job_problem(job_and_task, problem_state):
    if problem_state.start_time is not None:
        start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(problem_state.start_time))
    else:
        start_time = "start time not captured"
    return "%s\t%s\t%s.%s:\t%s" % \
           (start_time, problem_state.user, job_and_task.job_id, job_and_task.task_id, problem_state.state)

def format_problems(dct_problems):
    message = ""
    for host, problems in dct_problems.iteritems():
        message += "%s:\n\t%s\n" % (host, "\n\t".join(problems))
    return message

def sigint_handler(signum, frame):
    """Log termination cleanly"""
    logging.log_message("Exiting after received signal %d\n%s" %
                        (signum, "".join(traceback.format_list(traceback.extract_stack(frame)))))
    sys.exit()

def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--email", "-e", action="append",
                        help="email address(es) to which report should be sent.  Default: write report to stdout")
    parser.add_argument("--user", "-u", action="append",
                        help="User(s) whose jobs should be checked.  Default: all users")
    parser.add_argument("--user-file", type=read_list_file,
                        help="File containing list of users whose jobs should be checked, one per line.")
    parser.add_argument("--exclude-user", action="append",
                        help="User(s) whose jobs should not be checked.")
    parser.add_argument("--exclude-user-file", type=read_list_file,
                        help="File containing list of users whose jobs should not be checked, one per line.")
    parser.add_argument("--ignore-job", "-j", action="append",
                        help="Job ID(s) that should be ignored")
    parser.add_argument("--ignore-job-file", type=read_list_file,
                        help="File containing list of job ID(s) that should be ignored, one per line")
    parser.add_argument("--ignore-host", action="append",
                        help="Host(s) that should be ignored")
    parser.add_argument("--ignore-host-file", type=read_list_file,
                        help="File containing list of host(s) that should be ignored, one per line")
    parser.add_argument("--sleep-mins", "-s", type=float, default=60,
                        help="Time (in minutes) to sleep between status checks.  Default: %(default)s")
    parser.add_argument("--heartbeat-mins", type=float, default=24 * 60,
                        help="Time (in minutes) to sleep between sending message even if no problems.  "
                             "Default: %(default)s")
    parser.add_argument("--delete-mins", "-d", type=float, default=30,
                        help="A job in deletion state for more than this time (in minutes) is treated as a problem.  "
                             "Default: %(default)s")
    parser.add_argument("--transfer-mins", "-t", type=float, default=30,
                        help="A job in transfer state for more than this time (in minutes) is treated as a problem.  "
                             "Default: %(default)s")
    parser.add_argument("--cpu_fraction", "-c", type=float, default=0.0001,
                        help="A job with cpu time/wallclock time < this value is considered cpu-starved.  "
                             "Default: %(default)s")
    parser.add_argument("--starved-jobs-per-host", type=int, default=2,
                        help="A cpu-starved job is not reported unless there are at least this many cpu-starved jobs"
                             " on the host.  Multiple task for the same job are counted as one.  Default: %(default)s")
    parser.add_argument("--starved-task-array-hosts", type=int, default=4,
                        help="If a task array has cpu-starved jobs on at least this many hosts, it is assumed that "
                             "the job is waiting, rather than that the host is starving it of cpu.  Default: %(default)s")
    parser.add_argument("--wallclock_threshold", "-w", type=float, default=3,
                        help="A job is not considered cpu-starved until its wallclock time (minutes) is >= this value. "
                             "Default: %(default)s")
    parser.add_argument("--cpu-hog", type=float, default=2.0,
                        help="A job with (cpu time/wallclock time) - slots > this value is considered a cpu hog"
                             "Default: %(default)s")
    parser.add_argument("--verbose", "-v", action="count", default=0,
                        help="Write progress to stderr. Use this option multiple times for more verbosity")
    parser.add_argument("--smtp", default="smtp.broadinstitute.org",
                        help="SMTP host for sending mail. Default: %(default)s")
    parser.add_argument("--name", "-n", default="UGER",
                        help="Name of batch system to be used in messages.  Default: %(default)s")

    # key: host; value: time host failed ping
    dct_no_ping_hosts = {}

    # key: host; value: time bad state first seen
    dct_unreachable_queue_hosts = {}

    # key: JobAndTask; value: ProblemState(state, time state first seen, time reported)
    dct_strange_state_jobs = {}
    last_summary_time = time.time()

    signal.signal(signal.SIGINT, sigint_handler)

    options = parser.parse_args(args)
    # Do not send more than one exception message every 30 minutes, to avoid email avalanche
    mail_throttler = email_notify.MailThrottler(smtp=options.smtp, recipients=options.email, throttle=30 * 60,
                                                subject="UGER monitor unhandled exception")


    logging.log_message("Starting daemon")
    while True:
        try:
            # parse options again every time so that any files are reloaded
            options = parser.parse_args(args)
            logging.log_message("Options: " + str(options), 1, options.verbose)
            options.user = combine_option_lists(options.user, options.user_file, ["*"])
            excluded_users = frozenset(combine_option_lists(options.exclude_user, options.exclude_user_file))
            ignored_jobs = frozenset(combine_option_lists(options.ignore_job, options.ignore_job_file))
            ignored_hosts = frozenset(combine_option_lists(options.ignore_host, options.ignore_host_file))

            logging.log_message("Getting job list", 1, options.verbose)
            # remove jobs to be ignored
            jobs = [job for job in uge_functions.get_job_list(options.user)
                    if job.user not in excluded_users and
                    job.host is not None and
                    job.host not in ignored_hosts and job.job_and_task.job_id not in ignored_jobs]

            logging.log_message("Clearing strange jobs that are gone", 1, options.verbose)
            # Remove strange jobs that were not returned in job list
            current_job_set = set([job.job_and_task for job in jobs])
            for strange_job in dct_strange_state_jobs.keys():
                if strange_job not in current_job_set:
                    remove_strange_job(dct_strange_state_jobs, strange_job)

            # key: host; value: list of problem strings
            dct_problems = {}

            logging.log_message("Pinging hosts", 1, options.verbose)
            for newly_failed_host in ping_hosts(set(dct_no_ping_hosts.keys() + [job.host for job in jobs]),
                                                dct_no_ping_hosts):
                dct_problems.setdefault(newly_failed_host, []).append("did not respond to ping")

            logging.log_message("Checking queue states", 1, options.verbose)
            for newly_failed_host in check_host_queue_state(set(dct_unreachable_queue_hosts.keys() +
                                                                        [job.host for job in jobs]),
                                                            dct_unreachable_queue_hosts):
                dct_problems.setdefault(newly_failed_host, []).append("has queue in unreachable state")

            running_jobs = []
            logging.log_message("Checking job state", 1, options.verbose)
            transfer_state = "transfer"
            deletion_state = "deletion"
            for job in jobs:
                if is_transfer_state(job.state):
                    update_strange_jobs(dct_strange_state_jobs, job, transfer_state)
                elif is_delete_state(job.state):
                    update_strange_jobs(dct_strange_state_jobs, job, deletion_state)
                elif is_running_state(job.state):
                    running_jobs.append(job)
                    problem_state = dct_strange_state_jobs.get(job.job_and_task)
                    if problem_state is not None and \
                        (problem_state.state == transfer_state or problem_state.state == deletion_state):
                        logging.log_message("Job %s leaving strange state %s and now running" %
                                            (str(job), dct_strange_state_jobs[job.job_and_task]),
                                            1, options.verbose)
                        remove_strange_job(dct_strange_state_jobs, job.job_and_task)
                elif job.job_and_task in dct_strange_state_jobs:
                    # job not in running state anymore
                    logging.log_message("Odd job state " + str(job), 1, options.verbose)
                    remove_strange_job(dct_strange_state_jobs, job.job_and_task)

            logging.log_message("Checking cpu usage", 1, options.verbose)
            # Group by requires input is sorted in order to be grouped
            running_jobs.sort(cmp=lambda job1, job2: cmp((job1.job_and_task.job_id, job1.job_and_task.task_id),
                                                         (job2.job_and_task.job_id, job2.job_and_task.task_id)))
            for ignore_key, job_group in itertools.groupby(running_jobs, key=lambda j: j.job_and_task.job_id):
                # convert from iterable to list, because the first element is accessed twice
                jobs_in_group = list(job_group)
                if jobs_in_group[0].queue == "interactive":
                    # Interactive jobs can appear cpu-starved
                    continue
                job_id = jobs_in_group[0].job_and_task.job_id
                logging.log_message("Getting details for " + str(job_id), 2, options.verbose)
                dct_job_details = uge_functions.get_job_detail(job_id)
                if dct_job_details is None:
                    # all task for job appear to be gone
                    for job in jobs_in_group:
                        remove_strange_job(dct_strange_state_jobs, job.job_and_task)
                else:
                    for job in jobs_in_group:
                        job_details = dct_job_details.get(job.job_and_task.task_id)
                        if job_details is None:
                            # this job and task is done
                            remove_strange_job(dct_strange_state_jobs, job.job_and_task)
                        else:
                            cpu_fraction = cpu_starved(job_details, options.wallclock_threshold, options.cpu_fraction)
                            if cpu_fraction is not None:
                                update_strange_jobs(dct_strange_state_jobs, job, "cpu-starved %f" % cpu_fraction, job_details)
                            else:
                                cpu_hog_amount = cpu_hog(job_details, options.cpu_hog)
                                if cpu_hog_amount is not None:
                                    update_strange_jobs(dct_strange_state_jobs, job, "cpu-hog %f" % cpu_hog_amount, job_details)

            report_time_threshold = time.time() - 60 * options.sleep_mins
            deletion_time_threshold = time.time() - 60 * options.delete_mins
            transfer_time_threshold = time.time() - 60 * options.transfer_mins
            dct_count_starved_jobs_by_host, dct_starved_job_host_count = \
                count_starved_jobs_by_host(dct_strange_state_jobs,
                                           options.starved_jobs_per_host if options.verbose > 1 else None,
                                           options.starved_task_array_hosts if options.verbose > 1 else None)

            for job_and_task, problem_state in dct_strange_state_jobs.items():
                if problem_state.reported_time is None:
                    if problem_state.state == deletion_state:
                        report_problem = problem_state.first_seen_time <= deletion_time_threshold
                    elif problem_state.state == transfer_state:
                        report_problem = problem_state.first_seen_time <= transfer_time_threshold
                    elif problem_state.state.startswith("cpu-starved"):
                        report_problem = problem_state.first_seen_time <= report_time_threshold and \
                            dct_count_starved_jobs_by_host[problem_state.host] >= options.starved_jobs_per_host and \
                            dct_starved_job_host_count[job_and_task.job_id] < options.starved_task_array_hosts
                    else:
                        report_problem = problem_state.first_seen_time <= report_time_threshold
                    if report_problem:
                        dct_problems.setdefault(problem_state.host, []).append(format_job_problem(job_and_task, problem_state))
                        dct_strange_state_jobs[job_and_task] = problem_state._replace(reported_time=time.time())

            # Generate message for all problems even when only reporting new problems.
            do_heartbeat = time.time() - last_summary_time >= 60 * options.heartbeat_mins
            if len(dct_problems) > 0 or do_heartbeat:
                last_summary_time = time.time()
                all_problems = {}
                for host in dct_no_ping_hosts.keys():
                    all_problems.setdefault(host, []).append("did not respond to ping")
                for host in dct_unreachable_queue_hosts.keys():
                    all_problems.setdefault(host, []).append("has queue in unreachable state")
                for job_and_task, problem_state in dct_strange_state_jobs.items():
                    if problem_state.reported_time is not None:
                        all_problems.setdefault(problem_state.host, []). \
                            append(format_job_problem(job_and_task, problem_state))


            if len(dct_problems) > 0:
                message = "The following problems have been detected since the last check:\n\n" + \
                          format_problems(dct_problems) + \
                          "\n\nAll outstanding problems:\n\n" + format_problems(all_problems)
                logging.log_message(message)
                if options.email is not None:
                    email_notify.email_report(message, "New " + options.name + " problems", options.email, options.smtp)

            if do_heartbeat:
                message = "Previously reported problems that have not cleared:\n\n" + format_problems(all_problems)
                logging.log_message(message)
                if options.email is not None:
                    email_notify.email_report(message, "All " + options.name + " problems", options.email, options.smtp)

            logging.log_message("Sleeping for %s minutes" % options.sleep_mins, 1, options.verbose)
            time.sleep(options.sleep_mins * 60)
        except Exception, e:
            exc_str = traceback.format_exc()
            logging.log_message("Looping after unhandled exception: " + exc_str)
            mail_throttler.add_message(exc_str)


if __name__ == "__main__":
    sys.exit(main())
