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
Display statistics about the load on Grid Engine cluster: cores & RAM available, disabled, in use, queued, etc. 
"""
import argparse
import collections
import itertools
import sys
import logging
import tabularize
import time
import uge_functions
import histogram_tools

default_h_vmem = '1g'
default_queue = "broad"
JobResourceInfo = collections.namedtuple("JobResourceInfo",
                                         ["job_and_task", "project", "slots", "queue", "memory", "user", "host",
                                          "state", "submit_time", "start_time", "priority", "tickets", "h_rt"])


def parse_h_rt(h_rt):
    if h_rt is None:
        return None
    elif h_rt == 'infinity' or h_rt == 'INFINITY':
        return sys.maxint
    else:
        return int(h_rt)


def parse_running_job_resource_info(job_xml_tree):
    """
    Get info about all running jobs, including resource info (-r) and project (-ext)
    :param job_xml_tree: output of qstat -xml -j ...
    :return: list of JobResourceInfo
    """
    ancestor_finder = uge_functions.XmlAncestorFinder(job_xml_tree)
    lst_ret = []
    for task_number in job_xml_tree.findall("*//JAT_task_number"):
        task = ancestor_finder.get_ancestor(task_number)
        job = ancestor_finder.get_ancestor(task, 2)
        job_id = job.findtext('JB_job_number')
        task_id = task_number.text
        if task.find("JAT_granted_destin_identifier_list") is None:
            logging.log_message("Strange XML for running job (%s) task (%s)" % (job_id, task_id))
            continue
        slots_str = task.findtext('*//JG_slots')
        slots = int(slots_str)
        h_vmem_node = job.find("JB_hard_resource_list/element/[CE_name='h_vmem']")
        if h_vmem_node is None:
            memory = default_h_vmem
        else:
            memory = h_vmem_node.findtext("CE_stringval")
        h_rt_node = job.find("JB_hard_resource_list/element/[CE_name='h_rt']")
        h_rt = parse_h_rt(h_rt_node.findtext("CE_stringval")) if h_rt_node is not None else None
        queue = job.findtext('*//QR_name', default_queue)
        lst_ret.append(JobResourceInfo(job_and_task=uge_functions.JobAndTask(job_id=job_id, task_id=task_id),
                                       project=job.findtext('JB_project'),
                                       slots=slots,
                                       queue=queue,
                                       memory=slots * uge_functions.dehumanize_memory(memory),
                                       user=job.findtext('JB_owner'),
                                       host=task.findtext("*//JG_qhostname"),
                                       state="running",
                                       submit_time=float(job.findtext('JB_submission_time')) / 1000,
                                       start_time=float(task.findtext('JAT_start_time')) / 1000,
                                       priority=None,
                                       tickets=None,
                                       h_rt=h_rt))
    return lst_ret


def get_pending_job_resource_info(job_id, tree):
    """
    :param job_id: a single job ID to get resource info about
    :param tree: elementtree representation of job list.
    :return: a job resource info object.  Only one will be returned even if there are multiple tasks for the job
    """
    job = tree.find("*//job_list/[JB_job_number='%s']" % job_id)
    return parse_pending_job_resource_info(job)


def get_pending_jobs_resource_info(tree):
    return [parse_pending_job_resource_info(job) for job in tree.findall("*//job_list[@state='pending']")]


def parse_pending_job_resource_info(job):
    slots_str = job.find('slots')
    if slots_str is None:
        slots = 1
    else:
        slots = int(slots_str.text)
    queue = job.findtext('hard_req_queue', default_queue)
    memory = slots * uge_functions.dehumanize_memory(job.findtext("hard_request[@name='h_vmem']", default_h_vmem))
    h_rt_str = job.findtext("hard_request[@name='h_rt']")
    h_rt = parse_h_rt(h_rt_str)
    time_str = job.find('JB_submission_time').text
    # Strip trailing fractional seconds
    time_str = time_str[:time_str.find('.')]
    submit_time = time.mktime(time.strptime(time_str, '%Y-%m-%dT%H:%M:%S'))
    return JobResourceInfo(job_and_task=uge_functions.JobAndTask(job_id=job.find('JB_job_number').text, task_id=None),
                           project=job.findtext('JB_project'),
                           slots=slots,
                           queue=queue,
                           memory=memory,
                           user=job.findtext('JB_owner'),
                           host=None,
                           state=job.get('state'),
                           submit_time=submit_time,
                           start_time=None,
                           priority=float(job.findtext('JAT_prio')),
                           tickets=int(job.findtext('tickets')),
                           h_rt=h_rt)


def replace_tickets(jobs, tree):
    def get_tickets_for_job(job_id):
        job_node = tree.find("*//job_list/[JB_job_number='%s']" % job_id)
        tickets = job_node.findtext('tickets')
        return int(tickets) if tickets is not None else None

    return [job._replace(tickets=get_tickets_for_job(job.job_and_task.job_id)) for job in jobs]


ReportLine = collections.namedtuple("ReportLine",
                                    ["label", "slots", "memory", "shares", "pct_slots", "pct_memory", "pct_shares"])


def filter_none(lst): return [elem for elem in lst if elem is not None]


def make_report_line(label, hosts):
    hosts = list(hosts)
    return ReportLine(label=label,
                      slots=sum(filter_none([host.slots for host in hosts])),
                      memory=sum(filter_none([host.memory for host in hosts])),
                      shares=None,
                      pct_slots=None,
                      pct_memory=None,
                      pct_shares=None)


def group_and_report(iterable, key_fun, key_prefix="", key_filter=None):
    """
    Group the input elements according to the key, and generate a ReportLine for each group
    :param iterable: elements to be grouped and reported
    :param key_fun: function for extracting key from element
    :param key_prefix: prepended to the key in order to make ReportList.label
    :param key_filter: if not None, skip keys for which key_filter(key) == False
    :return: list of ReportLine
    """
    return [make_report_line(key_prefix + key, group) for key, group in
            itertools.groupby(sorted(iterable, key=key_fun), key_fun)
            if key_filter is None or key_filter(key)]


def report_by_queue_project_and_user(lst_jobs, label, dct_shares_by_project, top_n_memory_users=sys.maxint,
                                     top_n_slot_users=sys.maxint):
    lst_ret = [make_report_line("All " + label, lst_jobs)]
    lst_ret.extend(group_and_report(lst_jobs, key_fun=lambda job: job.queue, key_prefix="queue "))

    lst_by_project = [
        report_line._replace(label="project " + report_line.label, shares=dct_shares_by_project[report_line.label])
        for report_line in group_and_report(lst_jobs, key_fun=lambda job: job.project)]
    lst_ret.extend(lst_by_project)

    # Summarize by user
    user_summaries = group_and_report(lst_jobs, key_fun=lambda job: job.user, key_prefix="user ")
    # Sort users in descending order of memory use and slot use
    users_memory_sorted = sorted(user_summaries, key=lambda reportLine: reportLine.memory, reverse=True)
    users_slot_sorted = sorted(user_summaries, key=lambda reportLine: reportLine.slots, reverse=True)
    # report on the top users of memory or slots, sorted descending by slots
    users_sorted = sorted(frozenset(users_memory_sorted[:top_n_memory_users] + users_slot_sorted[:top_n_slot_users]),
                          key=lambda reportLine: reportLine.slots, reverse=True)
    lst_ret.extend(users_sorted)
    return lst_ret


time_format = "%m/%d/%Y %H:%M:%S"


def format_time(secs):
    return time.strftime(time_format, time.localtime(secs))


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--verbose", "-v", action="count", default=0,
                        help="Write progress to stderr. Use this option multiple times for more verbosity")
    parser.add_argument("--job", "-j", action="append",
                        help="Job ID(s) for which insight into pending state is desired")
    parser.add_argument("--pending-detail", action="store_true", default=False,
                        help="Print detail about jobs that have jumped in from of query job")
    parser.add_argument("--ignore-job", "-i", action="append", help="Job ID(s) to ignore when -j is specified.  "
                                                                    "This is sometimes necessary to work around bad UGE XML for some jobs",
                        default=[])
    parser.add_argument("--ignore-project", "-p", action="append",
                        help="Projects to ignore when finding all extant shares.  "
                             "This is sometimes necessary because of unused projects",
                        default=[])
    parser.add_argument("--ignore-queue", "-q", action="append", help="Queues to ignore when counting healthy hosts.  "
                                                                      "This is sometimes necessary because of unused queues",
                        default=[])
    parser.add_argument("--memory-users", "-m", type=int, default=10,
                        help="Print user stats for the top N users of memory.  Default: %(default)s")
    parser.add_argument("--slot-users", "-s", type=int, default=10,
                        help="Print user stats for the top N users of slots.  Default: %(default)s")
    parser.add_argument("--mem-histogram-binsize", type=float, default=10 * 1024 * 1024 * 1024.0,
                        help="Bin size for available memory histogram.  Default: %(default)s")
    parser.add_argument("--nocolor", default=False, action='store_true',
                        help="Do not use ANSI colors in output")
    parser.add_argument("--write-snapshot", "-w",
                        help="Write a tar file containing the output of the various Grid Engine "
                             "commands executed")
    parser.add_argument("--read-snapshot", "-r",
                        help="Instead of executing Grid Engine commands to obtain current state, "
                             "read snapshot tar file containing a past state")
    options = parser.parse_args(args)

    if options.read_snapshot is not None:
        logging.log_message("Loading Grid Engine snapshot from " + options.read_snapshot, level=1,
                            verbosity=options.verbose)
        uge_state = uge_functions.uge_snapshot(snapshot_tar=options.read_snapshot)
    else:
        logging.log_message("Taking snapshot of current Grid Engine state", level=1, verbosity=options.verbose)
        uge_state = uge_functions.uge_snapshot()

    if options.write_snapshot is not None:
        logging.log_message("Writing Grid Engine snapshot to " + options.write_snapshot, level=1,
                            verbosity=options.verbose)
        uge_state.save(options.write_snapshot)

    logging.log_message("Getting host info", level=1, verbosity=options.verbose)
    dct_hosts = uge_functions.parse_host_info(qhost_tree=uge_state.get_xml(uge_state.uge_snapshot_hosts),
                                              queues_tree=uge_state.get_xml(uge_state.uge_snapshot_queues),
                                              queues_to_ignore=options.ignore_queue)
    lst_hosts = dct_hosts.values()

    logging.log_message("Getting job info", level=1, verbosity=options.verbose)
    running_jobs = parse_running_job_resource_info(uge_state.get_xml(uge_state.uge_snapshot_running_jobs))

    all_hosts_report_line = make_report_line("all hosts", lst_hosts)
    lst_report = [all_hosts_report_line]
    lst_report.extend(group_and_report(lst_hosts, key_fun=lambda _host: _host.state, key_prefix="hosts in state ",
                                       key_filter=lambda state: state is not None and len(state) > 0))

    logging.log_message("Getting project info", level=1, verbosity=options.verbose)
    projects = [uge_state.get_project_info(project) for project in uge_state.get_projects() if
                project not in options.ignore_project]
    total_shares = sum([int(project.fshare) for project in projects])
    dct_shares_by_project = dict([(project.name, project.fshare) for project in projects])

    lst_report.append(make_report_line(
        "available resources on healthy hosts",
        [host._replace(slots=host.slots - host.slots_used, memory=host.memory - host.mem_used)
         for host in lst_hosts if host.state is None and host.slots_used is not None and host.mem_used is not None]))

    lst_report.extend(report_by_queue_project_and_user(running_jobs, "running jobs", dct_shares_by_project,
                                                       top_n_memory_users=options.memory_users,
                                                       top_n_slot_users=options.slot_users))

    # Add the percentages
    lst_report = [report._replace(pct_slots=int(report.slots * 100 / all_hosts_report_line.slots),
                                  pct_memory=int(report.memory * 100 / all_hosts_report_line.memory),
                                  pct_shares=None if report.shares is None else int(
                                      int(report.shares) * 100 / total_shares))
                  for report in lst_report]
    # Figure out which project(s) have gotten more slots or RAM than their percentage of shares
    # Leading Falses corresponds to the 2 header lines
    lst_hog = [False, False] + \
              [report.pct_shares is not None and (
                      report.pct_memory > report.pct_shares or report.pct_slots > report.pct_shares) for report in
               lst_report]

    print "# RESOURCES ALLOCATED TO RUNNING JOBS\n"
    lst_rows = [('class', 'slots', 'memory', 'pct_slots', 'pct_memory', 'pct_shares')] + \
               [(report.label, str(report.slots), uge_functions.humanize_memory(report.memory),
                 str(report.pct_slots) + '%',
                 str(report.pct_memory) + '%',
                 '' if report.pct_shares is None else str(report.pct_shares) + '%') for report in lst_report]

    if options.nocolor:
        clear_ansi_color = ''
        start_ansi_red = ''
    else:
        clear_ansi_color = '\033[0m'
        start_ansi_red = '\033[31m'

    for line, do_color in zip(tabularize.format_table(lst_rows, right_justify=[False, True, True, True, True, True]),
                              lst_hog):
        if do_color:
            print start_ansi_red + line + clear_ansi_color
        else:
            print line

    print "\n# SHARES BY PROJECT\n"
    lst_rows = [('project', 'shares', 'pct_shares')] + \
               [(project.name, str(project.fshare), str(int(int(project.fshare) * 100 / total_shares)) + '%')
                for project in projects]
    for line in tabularize.format_table(lst_rows, right_justify=[False, True, True]):
        print line

    print "\n# HISTOGRAM OF AVAILABLE SLOTS\n"
    slots_histogram = histogram_tools.make_histogram(
        [host.slots - host.slots_used for host in lst_hosts if host.slots_used is not None and host.state is None], 1)
    lst_rows = [("available_slots", "num_hosts", "num_hosts >= this bin")] + \
               [(str(slots), str(hosts), str(hosts_ge)) for slots, hosts, hosts_ge in slots_histogram]
    for line in tabularize.format_table(lst_rows, right_justify=[True, True, True], column_sep="   "):
        print line

    print "\n# HISTOGRAM OF AVAILABLE MEMORY\n"
    slots_histogram = histogram_tools.make_histogram(
        [host.memory - host.mem_used for host in lst_hosts if host.mem_used is not None and host.state is None],
        options.mem_histogram_binsize)
    lst_rows = [("available_memory", "num_hosts", "num_hosts >= this bin")] + \
               [(uge_functions.humanize_memory(mem, precision=0) + " <= mem < " +
                 uge_functions.humanize_memory(mem + options.mem_histogram_binsize, precision=0),
                 str(hosts), str(hosts_ge))
                for mem, hosts, hosts_ge in slots_histogram]
    for line in tabularize.format_table(lst_rows, right_justify=[True, True, True], column_sep="   "):
        print line

    if options.job is not None and len(options.job) > 0:
        logging.log_message("Getting start times of running jobs", level=1, verbosity=options.verbose)
        pending_jobs_resource_info = get_pending_jobs_resource_info(uge_state.get_xml(uge_state.uge_snapshot_jobs))
        for pending_job in options.job:
            resource_info_list = filter(lambda _job: _job.job_and_task.job_id == pending_job, pending_jobs_resource_info)
            if len(resource_info_list) == 0:
                print "\n# %s is not a pending job" % pending_job
                continue
            resource_info = resource_info_list[0]
            pending_job_submit_time = resource_info.submit_time
            queue_jumper_resource_infos = [job_info for job_info in running_jobs if
                                           job_info.submit_time > pending_job_submit_time]
            lst_report = report_by_queue_project_and_user(queue_jumper_resource_infos, "queue jumpers",
                                                          dct_shares_by_project,
                                                          top_n_memory_users=options.memory_users,
                                                          top_n_slot_users=options.slot_users)
            print "\n# QUEUE JUMPERS for job %s. user %s, queue %s, project %s, slots %d, memory %s, tickets %s, h_rt %s\n" % (
                resource_info.job_and_task.job_id, resource_info.user, resource_info.queue, resource_info.project,
                resource_info.slots,
                uge_functions.humanize_memory(resource_info.memory), resource_info.tickets,
                uge_functions.humanize_seconds(resource_info.h_rt)
            )
            lst_rows = [('class', 'slots', 'memory', 'pct_slots', 'pct_memory')] + \
                       [(report.label, str(report.slots), uge_functions.humanize_memory(report.memory),
                         str(int(report.slots * 100 / all_hosts_report_line.slots)) + '%',
                         str(int(report.memory * 100 / all_hosts_report_line.memory)) + '%') for report in lst_report]
            for line in tabularize.format_table(lst_rows, right_justify=[False, True, True, True, True]):
                print line

            print "\n# HISTOGRAM OF SLOT RESERVATIONS FOR QUEUE JUMPERS for job %s\n" % pending_job
            queue_jumper_slots_histo = histogram_tools.make_histogram(
                [job_info.slots for job_info in queue_jumper_resource_infos], 1)
            lst_rows = [("slots_used", "num_queue_jumper_jobs", "num_queue_jumper_jobs >= this bin")] + \
                       [(str(slots), str(jobs), str(jobs_ge)) for slots, jobs, jobs_ge in queue_jumper_slots_histo]
            for line in tabularize.format_table(lst_rows, right_justify=[True, True, True], column_sep="   "):
                print line

            lst_available_hosts = [(host.host, host.slots - host.slots_used, host.memory - host.mem_used)
                                   for host in lst_hosts
                                   if host.state is None and host.mem_used is not None and
                                   host.slots - host.slots_used >= resource_info.slots and
                                   host.memory - host.mem_used >= resource_info.memory and
                                   resource_info.queue in host.queues]
            # Sort by (available_slots, available_memory)
            lst_available_hosts.sort(key=lambda _tup: _tup[1:])
            # Convert to strings for printing
            lst_available_hosts = [(tup[0], str(tup[1]), uge_functions.humanize_memory(tup[2])) for tup in
                                   lst_available_hosts]
            print "\n# %d HEALTHY HOSTS WITH AT LEAST %d SLOTS AND %s MEMORY IN %s QUEUE FOR JOB %s\n" % \
                  (len(lst_available_hosts), resource_info.slots, uge_functions.humanize_memory(resource_info.memory),
                   resource_info.queue, pending_job)
            lst_rows = [("host", "available_slots", "available_memory")] + lst_available_hosts
            for line in tabularize.format_table(lst_rows, right_justify=[False, True, True], column_sep="   "):
                print line

            # report on pending jobs with more tickets than query job
            lst_high_ticket_pending_jobs = [job_info for job_info in pending_jobs_resource_info if
                                            job_info.tickets > resource_info.tickets]

            print "\n# Pending jobs with more tickets than job " + pending_job
            lst_rows = [('class', 'slots', 'memory')] + \
                       [(report.label, str(report.slots), uge_functions.humanize_memory(report.memory)) for report in
                        (report_by_queue_project_and_user(lst_high_ticket_pending_jobs, "pending jobs",
                                                          dct_shares_by_project))]
            for line in tabularize.format_table(lst_rows, right_justify=[False, True, True]):
                print line

            print "\n# Pending jobs with more tickets than job " + pending_job + " but submitted after"
            lst_high_ticket_pending_jumpers = [job_info for job_info in lst_high_ticket_pending_jobs if
                                               job_info.submit_time > resource_info.submit_time]
            lst_rows = [('class', 'slots', 'memory')] + \
                       [(report.label, str(report.slots), uge_functions.humanize_memory(report.memory)) for report in
                        (report_by_queue_project_and_user(lst_high_ticket_pending_jumpers,
                                                          "pending high-ticket queue jumpers", dct_shares_by_project))]
            for line in tabularize.format_table(lst_rows, right_justify=[False, True, True]):
                print line

            if options.pending_detail:
                # get tickets for running jobs
                queue_jumper_resource_infos = replace_tickets(queue_jumper_resource_infos,
                                                              uge_state.get_xml(uge_state.uge_snapshot_jobs))
                # Both these tables are in descending order by tickets
                queue_jumper_resource_infos.sort(cmp=lambda x, y: cmp(y.tickets, x.tickets))
                lst_high_ticket_pending_jumpers.sort(cmp=lambda x, y: cmp(y.tickets, x.tickets))
                lst_header = ['job', 'task', 'user', 'project', 'slots', 'memory', 'tickets', 'h_rt', 'submitted',
                              'dispatched']
                print "\n# Running queue jumper details"
                lst_rows = [lst_header] + [(job.job_and_task.job_id, job.job_and_task.task_id, job.user, job.project,
                                            str(job.slots), uge_functions.humanize_memory(job.memory), str(job.tickets),
                                            uge_functions.humanize_seconds(job.h_rt),
                                            format_time(job.submit_time), format_time(job.start_time))
                                           for job in queue_jumper_resource_infos]
                for line in tabularize.format_table(lst_rows, right_justify=[False, True, False, False, True, True,
                                                                             True, True, False, False]):
                    print line

                lst_header = ['job', 'user', 'project', 'slots', 'memory', 'tickets', 'h_rt', 'submitted']
                print "\n# Pending queue jumper details"
                lst_rows = [lst_header] + [(job.job_and_task.job_id, job.user, job.project,
                                            str(job.slots), uge_functions.humanize_memory(job.memory), str(job.tickets),
                                            uge_functions.humanize_seconds(job.h_rt), format_time(job.submit_time))
                                           for job in lst_high_ticket_pending_jumpers]
                for line in tabularize.format_table(lst_rows, right_justify=[False, False, False, True, True,
                                                                             True, True, False]):
                    print line


if __name__ == "__main__":
    sys.exit(main())
