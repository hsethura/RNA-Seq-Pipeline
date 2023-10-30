#!/usr/bin/env python
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2017 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.
"""
Monitor Grid Engine accounting file for 'black holes,' i.e. hosts to which many jobs are scheduled and terminate
abnormally and quickly over a short period of time.  This script runs until killed.
"""
import argparse
import collections
import logging
import subprocess
import sys
import time

import os
import email_notify

SGE_ROOT = 'SGE_ROOT'
SGE_CELL = 'SGE_CELL'
accounting_subpath = 'common/accounting'
if SGE_ROOT in os.environ and SGE_CELL in os.environ:
    full_accounting_path = os.path.join(os.environ[SGE_ROOT], os.environ[SGE_CELL], accounting_subpath)
else:
    full_accounting_path = None

AccountingLine = collections.namedtuple("AccountingLine",
                                        ["job_id",          # string
                                         "task_number",     # int
                                         "host",            # string
                                         "owner",            # string
                                         "start_time",      # float (seconds)
                                         "end_time",        # float (seconds)
                                         "exit_status",     # int
                                         "failed",          # string
                                         "duration",        # float (seconds),
                                         "deleted_by"      # string
                                         ])

AccountingLine.__str__ = lambda self: "job=%s, task=%d, host=%s, owner=%s, start_time=%s, end_time=%s, " \
                                      "exit_status=%d, failed=%s, duration=%f, deleted_by=%s" % (
    self.job_id, self.task_number, self.host, self.owner, time.ctime(self.start_time), time.ctime(self.end_time),
    self.exit_status, self.failed, self.duration, self.deleted_by)


def parse_accounting_line(str_line, time_conversion_factor=1000.0):
    """
    :param str_line: colon-separated accounting line, with trailing newline removed
    :param time_conversion_factor: For SGE, pass 1
    :return: AccountingLine named tuple
    """
    lst_fields = str_line.split(':')
    # SGE has 45 fields, UGE has 53, but first 45 are the same as for SGE
    if len(lst_fields) < 45:
        raise Exception("Not enough fields (%d) in accounting line: %s" % (len(lst_fields), str_line))
    host = lst_fields[1]
    owner = lst_fields[3]
    job_id = lst_fields[5]
    task_number = lst_fields[35]
    if len(task_number) > 0:
        task_number = int(task_number)
    else:
        task_number = None
    failed = lst_fields[11]
    exit_status = int(lst_fields[12])
    # Convert times to seconds
    start_time = float(int(lst_fields[9]) / time_conversion_factor)
    end_time = float(int(lst_fields[10]) / time_conversion_factor)
    duration = end_time - start_time
    deleted_by = None
    if len(lst_fields) > 46:
        # accounting manpage is wrong about the position of this value
        deleted_by = lst_fields[46]
        if deleted_by == 'NONE':
            deleted_by = None
    return AccountingLine(job_id=job_id, task_number=task_number, host=host, start_time=start_time, end_time=end_time,
                          exit_status=exit_status, failed=failed, duration=duration, deleted_by=deleted_by, owner=owner)


def generate_accounting_lines(f_tail):
    while True:
        yield parse_accounting_line(f_tail.readline().rstrip("\n"))


def update_dict(d, new_job, age_out_interval):
    """
    Only keep most recent job for a given user on a host, in order to avoid thinking that a host is a black hole
    because a single user has submitted a lot of defective jobs.
    """
    if new_job.host not in d:
        d[new_job.host] = {new_job.owner: new_job}
    else:
        dct_jobs_for_host = d[new_job.host]
        for k, job in dct_jobs_for_host.items():
            if new_job.end_time - job.end_time > age_out_interval:
                del dct_jobs_for_host[k]
        dct_jobs_for_host[new_job.owner] = new_job


def format_failed_jobs_dict(dct_by_host):
    return "\n".join([host + "\n   " + "\n   ".join([str(job) for job in jobs_for_host.itervalues()])for host, jobs_for_host in dct_by_host.iteritems()])

def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-duration", "-j", type=float, default=60.0,
                        help="Count jobs that die within this many seconds.  Default: %(default)s")
    parser.add_argument("--num-failed-jobs", "-n", type=int, default=5,
                        help="Report if at least this many distinct jobs "
                             "(as opposed to multiple tasks for the same job)"
                             " die on the same host.  Default: %(default)s")
    parser.add_argument("--interval", "-i", type=float, default=300.0,
                        help="Report if enough jobs failed on the same host within this many seconds.  "
                             "Default: %(default)s")
    if full_accounting_path is None:
        parser.add_argument("--accounting", "-a", required=True,
                            help="Grid Engine accounting file to monitor.  "
                                 "This is required if SGE_ROOT and SGE_CELL environment variables are not set.  "
                                 "Load a Grid Engine dotkit to set these environment variables.")
    else:
        parser.add_argument("--accounting", "-a", default=full_accounting_path,
                            help="Grid Engine accounting file to monitor.  Default: %(default)s")

    parser.add_argument("--email", "-e", action="append",
                        help="email address(es) to which report should be sent.  Default: only write report to stdout")
    parser.add_argument("--smtp", default="smtp.broadinstitute.org",
                        help="SMTP host for sending mail. Default: %(default)s")
    parser.add_argument("--email-interval", type=int, default=300,
                        help="Send an email message no more frequently than this many seconds, in order to avoid "
                             "swamping mail server.  Default: %(default)s")
    parser.add_argument("--heartbeat-interval", type=int, default=60*60*24,
                        help="Send an email message after this many seconds, to make sure the daemon is alive.  "
                             "Default: %(default)s")
    parser.add_argument("--verbose", "-v", action="count", default=0,
                        help="Write progress to stderr. Use this option multiple times for more verbosity")
    options = parser.parse_args(args)
    logging.log_message(str(options))

    mail_throttler = email_notify.MailThrottler(smtp=options.smtp, recipients=options.email, throttle=options.email_interval,
                                   subject="UGER black hole hosts")

    last_heartbeat = time.time()

    # bufsize=-1 => 'fully buffered'
    proc = subprocess.Popen(['tail', '-F', options.accounting], stdout=subprocess.PIPE, bufsize=-1)
    # Skip first line which may be truncated
    proc.stdout.readline()
    dct_by_host = dict()
    for job in generate_accounting_lines(proc.stdout):
        now = time.time()
        if now - last_heartbeat >= options.heartbeat_interval:
            email_notify.email_report(message="Grid Engine black hole detector is alive",
                                      subject="Grid Engine black hole heartbeat",
                                      recipients=options.email, smtp_host=options.smtp)
            logging.log_message("Grid Engine black hole detector is alive")
            last_heartbeat = now

        mail_throttler.maybe_send()
        logging.log_message("Terminated job: " + str(job), level=3, verbosity=options.verbose)
        if job.exit_status >= 128 and job.duration <= options.job_duration and job.deleted_by is None:
            b_potential_black_hole = False
            logging.log_message("Failed job " + str(job), level=2, verbosity=options.verbose)
            update_dict(dct_by_host, job, age_out_interval=options.interval * 2)
            for host in dct_by_host.iterkeys():
                dct_host = dct_by_host[host]
                if len(dct_host) >= options.num_failed_jobs:
                    b_potential_black_hole = True
                    logging.log_message("potential black hole " + host, level=1, verbosity=options.verbose)
                    sorted_jobs = sorted(dct_host.itervalues(), key=lambda j: j.end_time)
                    for i in xrange(0, len(sorted_jobs) - options.num_failed_jobs):
                        if sorted_jobs[i+options.num_failed_jobs].end_time - sorted_jobs[i].end_time <= \
                                options.interval:
                            black_hole_message = "BLACK HOLE: %s has %d failed jobs" % (host, len(sorted_jobs))
                            logging.log_message(black_hole_message, file=sys.stdout)
                            logging.log_message(black_hole_message)
                            mail_throttler.add_message(black_hole_message)
                            break
            if b_potential_black_hole:
                logging.log_message("Current failed job state: " + format_failed_jobs_dict(dct_by_host), level=2,
                                    verbosity=options.verbose)



if __name__ == "__main__":
    sys.exit(main())
