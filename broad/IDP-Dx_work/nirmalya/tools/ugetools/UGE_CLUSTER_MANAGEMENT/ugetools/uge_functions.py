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
Note that files containing lists of users, hosts, jobs are re-read every time status is checked, so these files
can be changed without having to restart the daemon.
"""
import collections
import itertools
import re
import subprocess
import time
import xml.etree.ElementTree
import logging
import tarfile
import os
import StringIO
import curses.ascii

class ProcessException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

def run_command_capture_output(lst_args):
    """
    :param lst_args: Command to invoke plus args to pass 
    :return: command stdout as a string
    """
    proc = subprocess.Popen(lst_args, stdout=subprocess.PIPE)

    str_stdout = proc.stdout.read()
    if proc.wait() != 0:
        raise ProcessException("ERROR: " + lst_args[0] + " returned " + str(proc.returncode))
    return str_stdout


def run_command_parse_xml(lst_args):
    """
    :param lst_args: Command to invoke plus args to pass 
    :return: tuple(aElementTree object with the results, xml as a string)
    """
    str_xml = run_command_capture_output(lst_args)
    tree = parse_xml("Command: '%s'" % ' '.join(lst_args), str_xml)
    return tree, str_xml


# Backslash needs to be considered a bad character
bad_xml_char = re.compile("[^ -[\]-~\s]")
# These values cannot be converted into entity references
st_bad_entities = frozenset([27])

def clean_bad_chars(str_xml):
    """
    Replace invalid ASCII with entity references.  Grid Engine doesn't clean up non-printable ASCII, e.g. from
    environment variables.
    """
    matches = [mo for mo in bad_xml_char.finditer(str_xml)]
    matches.reverse()
    for match in matches:

        numeric_value = ord(match.group())
        replacement = "?" if numeric_value in st_bad_entities else "&#%d;" % numeric_value
        str_xml = str_xml[:match.start()] + replacement + str_xml[match.end():]
    return str_xml

def parse_xml(message, str_xml_raw):
    """
    Parse xml string and report errors nicely
    :param message: used only in exception strings
    :param str_xml: xml string to be parsed
    :return: ElementTree
    """
    try:
        str_xml = clean_bad_chars(str_xml_raw)
        tree = xml.etree.ElementTree.fromstring(str_xml)
    except xml.etree.ElementTree.ParseError, e:
        # When reporting bad XML, show 200 lines of context
        lines = str_xml.split('\n')
        bad_lineno = e.position[0] - 1
        context = lines[max(0, bad_lineno - 100): bad_lineno] + [">>>" + lines[bad_lineno]] + lines[
                                                                                              bad_lineno + 1:bad_lineno + 100]
        raise Exception("%r: %s\n%r" % (e, message, "\n".join(context)))
    except Exception, e:
        raise Exception("%r: %s\n%r" % (e, message, str_xml[:1000]))
    return tree


class XmlAncestorFinder(object):
    def __init__(self, tree):
        self.__dct_parent = {c:p for p in tree.iter() for c in p}

    def get_ancestor(self, node, degree=1):
        if degree < 0:
            raise Exception("ancestor degree(%s) cannot be < 0" % degree)
        elif degree == 0:
            return node
        else:
            return self.get_ancestor(self.__dct_parent[node], degree-1)


def make_xml_parent_map(tree):
    return {c:p for p in tree.iter() for c in p}

def run_qstat(lst_args):
    """
    :param lst_args: Invoke qstat with -xml plus these args
    :return: tuple(aElementTree object with the results, xml as a string)
    """
    return run_command_parse_xml(["qstat", "-xml"] + lst_args)


def run_qhost(lst_args):
    """
    :param lst_args: Invoke qhost with -xml plus these args
    :return: tuple(aElementTree object with the results, xml as a string)
    """
    return run_command_parse_xml(["qhost", "-xml"] + lst_args)


def get_host_from_queue_name(queue_name):
    fields = queue_name.split("@", 1)
    if len(fields) == 1:
        return queue_name
    else:
        return fields[1]

def get_queue_from_queue_name(queue_name):
    fields = queue_name.split("@", 1)
    if len(fields) == 1:
        return queue_name
    else:
        return fields[0]

# Both of these are strings.  task_id may be None for a non-running job
JobAndTask = collections.namedtuple("JobAndTask", ["job_id", "task_id"])

JobInfo = collections.namedtuple("JobInfo", ["job_and_task", "state", "user", "host", "queue"])

def parse_time(strTime):
    """Parse a string like this: 2017-05-11T10:02:03.589"""
    # Strip fractional seconds.
    strTime = strTime.split(".")[0]
    return time.mktime(time.strptime(strTime, '%Y-%m-%dT%H:%M:%S'))

def get_job_list(users, other_args=None):
    """
    :param users: get jobs for these users
    :param other_args: additional qstat args
    :return: List of JobInfo objects.  host and task_id may be None for non-running jobs
    """
    if other_args is None:
        other_args = []
    tree, str_xml = run_qstat([val for val in itertools.chain.from_iterable([['-u', val] for val in users])] +
                              other_args)

    running_jobs = [(job_listElement.find("JB_job_number").text,
                     job_listElement.find("state").text,
                     job_listElement.findtext("tasks"),
                     job_listElement.find("JB_owner").text,
                     job_listElement.findtext("queue_name")) for queue_infoElement
                    in tree.findall("queue_info")
                    for job_listElement in queue_infoElement.findall("job_list")]
    running_jobs = [JobInfo(JobAndTask(job_id, task_id or "1"), state, user,
                            get_host_from_queue_name(queue_name), get_queue_from_queue_name(queue_name))
                    for job_id, state, task_id, user, queue_name in running_jobs]

    non_running_jobs = [JobInfo(JobAndTask(job_listElement.find("JB_job_number").text, None),
                                job_listElement.find("state").text,
                                job_listElement.find("JB_owner").text,
                                None, None)
                        for job_infoElement in tree.findall("job_info")
                        for job_listElement in job_infoElement.findall("job_list")]

    return running_jobs + non_running_jobs


JobDetail = collections.namedtuple("JobDetail", ["job_and_task", "user", "job_name",
                                                 "submission_time_secs", "start_time_secs",
                                                 "cpu_secs", "slots"])


def int_or_none(element):
    if element is None:
        return None
    else:
        return int(element.text)


def get_job_detail(job_id):
    """
    :param job_id: Either a job_id string, or job_id.task_array_index string
    :return: dictionary key: task_id, value: JobDetail object for the running task, or None if job is gone
    """
    try:
        tree, str_xml = run_qstat(['-j', job_id])
    except Exception, e:
        logging.log_message("Exception getting job detail for " + job_id + ": " + str(e))
        return None
    if tree.tag == 'unknown_jobs':
        # Assume job has completed
        return None

    ret = parse_job_detail(tree.find('djob_info/element'))

    return ret


def parse_job_detail(djob_info):
    user = djob_info.find('JB_owner').text
    job_id = djob_info.find('JB_job_number').text
    job_name = djob_info.find('JB_job_name').text
    submission_time_secs = int(djob_info.find('JB_submission_time').text) / 1000.0
    slots_element = djob_info.find(".//JG_slots")
    if slots_element is None:
        slots = None
    else:
        slots = int(slots_element.text)
    ret = {}
    for element in djob_info.findall(".//JB_ja_tasks/element"):
        task_id = element.find("JAT_task_number").text
        start_time_secs = int_or_none(element.find('JAT_start_time'))  # None if not started
        if start_time_secs is not None:
            start_time_secs = start_time_secs / 1000.0
        for event in element.findall(".//JAT_scaled_usage_list/Events"):
            if event.find("UA_name").text == "cpu":
                cpu_secs = float(event.find("UA_value").text)
                break
        else:
            cpu_secs = 0.0

        ret[task_id] = JobDetail(job_and_task=JobAndTask(job_id=job_id, task_id=task_id),
                                 user=user,
                                 job_name=job_name,
                                 submission_time_secs=submission_time_secs,
                                 start_time_secs=start_time_secs,
                                 cpu_secs=cpu_secs,
                                 slots=slots)
    return ret

def get_multiple_running_job_detail(lst_job_id):
    """
    :param lst_job_id: a list of jobs to query 
    :return: dictionary key: job_id, value: dictionary with key: task_id, value: JobDetail
    """
    try:
        jobs = ",".join(lst_job_id)
        tree, str_xml = run_qstat(['-j', jobs])
    except Exception, e:
        logging.log_message("Exception getting job detail for " + jobs + ": " + str(e))
        return None
    if tree.tag == 'unknown_jobs':
        # Assume job has completed
        return None

    lst_dicts = [parse_job_detail(element) for element in tree.findall('djob_info/element')]

    return dict([(dct_job.itervalues().next().job_and_task.job_id, dct_job)for dct_job in lst_dicts])

def get_hosts_with_unreachable_queues(hosts):
    """
    Query all the queues on all the hosts (qstat -f -q \*@host)
    :param hosts: Hosts to be queried
    :return: set of hosts that has at least one queue with a bad state
    """
    args = ["-f"]
    for host in hosts:
        args += ["-q", "*@" + host]
    tree, str_xml = run_qstat(args)
    ret = set()
    for queue in tree.findall('./queue_info/Queue-List'):
        state_element = queue.find('state')
        if state_element is not None and state_element.text.find("u") != -1:
            ret.add(queue.find('name').text.split("@", 1)[1])

    return ret


HostInfo = collections.namedtuple("HostInfo", ["host", "slots", "memory", "state", "slots_used", "mem_used", "queues"])


def get_host_info():
    """
    :return: dictionary key: host, value HostInfo 
    """
    qhost_tree, qhost_str_xml = run_qhost(['-q'])
    qstat_tree, qstat_str_xml = run_qstat(['-f'])
    return parse_host_info(qhost_tree, qstat_tree)

# Each element is a set of states that are not considered in conflict with one another.
lst_compatible_host_states = [
    frozenset(['u', 'uE', 'd'])
]

def is_host_state_compatible(s1, s2):
    if s1 == s2: return True
    for st_compatible_states in lst_compatible_host_states:
        if s1 in st_compatible_states:
            return s2 in st_compatible_states
    return False


def parse_host_info(qhost_tree, queues_tree, queues_to_ignore=[]):
    """
    :return: dictionary key: host, value HostInfo
    """
    dctRet = {}
    for host_node in qhost_tree.findall('host'):
        host_name = host_node.get('name')
        dct_hostvalues = dict([(hostvalue_node.get('name'), hostvalue_node.text) for hostvalue_node in host_node.findall('hostvalue')])
        if dct_hostvalues['num_proc'] != '-':
            slots = int(dct_hostvalues['num_proc'])
            slots_used = sum([int(slots_used_node.text) for slots_used_node in host_node.findall(".//queuevalue[@name='slots_used']")])
            memory = dehumanize_memory(dct_hostvalues['mem_total'])
            mem_used = 0 if dct_hostvalues['mem_used'] == '-' else dehumanize_memory(dct_hostvalues['mem_used'])
            dctRet[host_name] = HostInfo(host=host_name, slots=slots, memory=memory, state=None, slots_used=slots_used,
                                         mem_used=mem_used, queues=set())
        else:
            dctRet[host_name] = HostInfo(host=host_name, slots=None, memory=None, state=None, slots_used=None,
                                         mem_used=None, queues=set())
    for queue_info in queues_tree.findall('*/Queue-List'):
        state = queue_info.findtext('state')
        if state is None: state = ''
        # Ignore suspended state
        state = re.sub('s', '', state)
        # Ignore configuration ambiguous state
        state = re.sub('c', '', state)
        # If disabled, ignore other state flags, because they can vary between queues on a host
        if 'd' in state:
            state = 'd'

        queue = queue_info.findtext('name')
        queue_split = queue.split('@', 1)
        host = queue_split[1]
        queue_name = queue_split[0]
        if queue_name in queues_to_ignore:
            continue
        host_info = dctRet.get(host)
        host_info.queues.add(queue_name)

        if len(state) > 0:
            if host_info is None:
                logging.log_message(host + " found in qstat but not qhost")
            elif host_info.state is None:
                dctRet[host] = host_info._replace(state=state)
            elif not is_host_state_compatible(host_info.state, state):
                raise Exception("Conflicting states for %s: %s != %s" % (host, host_info.state, state))

    return dctRet

dctNumericSuffix = {
    "k": 1024,
    "m": 1024 * 1024,
    "g": 1024 * 1024 * 1024,
    "": 1
}

# Create list for converting from bytes to number with suffix
lstHumanize = zip(dctNumericSuffix.values(), dctNumericSuffix.keys())
lstHumanize.sort(reverse=True)

class MemoryParseError(Exception):
    pass


def dehumanize_memory(strMem):
    """Returns an integer number of bytes, by translating suffices like 'm', 'k', 'g'"""
    if strMem == "-":
        return float(0)
    elif strMem[-1].isalpha():
        strSuffix = strMem[-1].lower()
        strNumber = strMem[:-1]
    else:
        strNumber = strMem
        strSuffix = ""

    if strSuffix not in dctNumericSuffix:
        raise MemoryParseError("Unrecognized memory suffix: " + strSuffix)
    try:
        return float(strNumber) * dctNumericSuffix[strSuffix]
    except ValueError, v:
        raise MemoryParseError("Problem parsing memory specification '%s': %s" % (strMem, v.message))

def humanize_memory(fMem, precision=1):
    final_suffix = " "
    for (divisor, suffix) in lstHumanize:
        if fMem > divisor:
            final_suffix = suffix
            fMem = fMem/float(divisor)
            break
    return ("%." + str(precision) + "f%s") % (fMem, suffix)

MinuteInSeconds = 60
HourInSeconds = 60 * MinuteInSeconds
DayInSeconds = 24 * HourInSeconds
def humanize_seconds(seconds):
    if seconds is None:
        return None
    seconds = int(seconds)
    days = seconds / DayInSeconds
    remainder = seconds % DayInSeconds
    hours = remainder/HourInSeconds
    remainder = remainder % HourInSeconds
    minutes = remainder/MinuteInSeconds
    remainder = remainder % MinuteInSeconds
    if days != 0:
        return "%d:%02d:%02d:%02d" % (days, hours, minutes, remainder)
    return "%d:%02d:%02d" % (hours, minutes, remainder)

def get_projects():
    str_out = run_command_capture_output(["qconf", "-sprjl"])
    return str_out.splitlines()

ProjectInfo = collections.namedtuple("ProjectInfo", ["name", "oticket", "fshare", "acl", "xacl"])

def get_project_info(project):
    return parse_project_info(run_command_capture_output(["qconf", "-sprj", project]))


def parse_project_info(str_out):
    lines = str_out.splitlines()
    dct_project = dict([tuple(line.split(' ')) for line in lines])
    return ProjectInfo(**dct_project)


class uge_snapshot(object):
    uge_snapshot_start = "start.txt"
    uge_snapshot_end = "end.txt"
    uge_snapshot_hosts = "hosts.xml"
    uge_snapshot_queues = "queues.xml"
    uge_snapshot_jobs = "jobs.xml"
    uge_snapshot_running_jobs = "running_jobs.xml"
    uge_snapshot_projects = "projects.txt"
    uge_snapshot_project_prefix = "project"
    uge_snapshot_command_suffix = ".cmd"

    def __get_subdict(self, lst_parent_keys):
        dct = self.dct
        for key in lst_parent_keys: dct = dct[key]
        return dct

    def __snapshot_one_command(self, output_label, args, lst_parent_keys = []):
        """Capture command string, and command output"""
        dct = self.__get_subdict(lst_parent_keys)
        dct[output_label] = run_command_capture_output(args)
        dct[output_label + uge_snapshot.uge_snapshot_command_suffix] = " ".join(args)

    def __init__(self, snapshot_tar=None, lst_users=None):
        self.dct = {}
        self.xml_cache = {}
        if snapshot_tar is None:
            self.dct[uge_snapshot.uge_snapshot_start]  = time.ctime()
            self.__snapshot_one_command(uge_snapshot.uge_snapshot_hosts, ["qhost", "-xml", "-q"])
            self.__snapshot_one_command(uge_snapshot.uge_snapshot_queues, ["qstat", "-xml", "-f"])
            users = lst_users if lst_users is not None else ['*']
            lst_user_args = [arg for user in users for arg in ["-u", user]]
            self.__snapshot_one_command(uge_snapshot.uge_snapshot_jobs, ["qstat", "-xml", "-ext", "-r", "-pri", "-urg"] + lst_user_args)
            self.__snapshot_one_command(uge_snapshot.uge_snapshot_projects, ["qconf", "-sprjl"])

            # this XML parsing tends not to fail, so don't worry about handling exceptions
            xml_job_list = self.get_xml(uge_snapshot.uge_snapshot_jobs)
            lst_running_jobs = [node.text for node in xml_job_list.findall("*/job_list[@state='running']/JB_job_number")]
            self.__snapshot_one_command(uge_snapshot.uge_snapshot_running_jobs, ["qstat", "-xml", "-j", ",".join(lst_running_jobs)])

            self.__snapshot_one_command(uge_snapshot.uge_snapshot_projects, ["qconf", "-sprjl"])
            self.dct[uge_snapshot.uge_snapshot_project_prefix] = {}
            for project in self.dct[uge_snapshot.uge_snapshot_projects].splitlines():
                self.__snapshot_one_command(project, ["qconf", "-sprj", project], lst_parent_keys=[uge_snapshot.uge_snapshot_project_prefix])

            self.dct[uge_snapshot.uge_snapshot_end]  = time.ctime()
        else:
            if lst_users is not None:
                raise Exception("Cannot specify users when loading snapshot")
            tf = tarfile.open(snapshot_tar)
            for ti in tf.getmembers():
                if ti.isdir():
                    continue
                content = tf.extractfile(ti).read()
                ti_path = ti.name.split("/")
                self.__load_uge_snapshot_helper(self.dct, ti_path, content)

            tf.close()

    def get_xml(self, label, lst_parent_keys=[]):
        # Cache parsed XML
        cache_key = list(lst_parent_keys)
        cache_key.append(label)
        cache_key = tuple(cache_key)
        if cache_key not in self.xml_cache:
            dct = self.__get_subdict(lst_parent_keys)
            self.xml_cache[cache_key] = parse_xml(dct[label + uge_snapshot.uge_snapshot_command_suffix], dct[label])
        return self.xml_cache[cache_key]

    def get(self, label, lst_parent_keys=[]):
        dct = self.__get_subdict(lst_parent_keys)
        return dct[label]

    def get_projects(self):
        return sorted(filter(lambda key: not key.endswith(uge_snapshot.uge_snapshot_command_suffix),
                      self.dct[uge_snapshot.uge_snapshot_project_prefix].keys()))

    def get_project_info(self, project):
        return parse_project_info(self.get(project, lst_parent_keys=[uge_snapshot.uge_snapshot_project_prefix]))

    def save(self, path):
        mode = "w"
        if path.endswith(".gz"):
            mode += ":gz"
        tf = tarfile.open(path, mode)
        self.__save_uge_snapshot_helper(tf, self.dct, time.time())
        tf.close()

    def __save_uge_snapshot_helper(self, tf, dct, mtime, prefix=""):
        for key, val in dct.iteritems():
            if isinstance(val, str):
                ti = tarfile.TarInfo(os.path.join(prefix, key))
                ti.size = len(val)
                ti.mtime = mtime
                tf.addfile(ti, StringIO.StringIO(val))
            elif isinstance(val, dict):
                self.__save_uge_snapshot_helper(tf, val, mtime, os.path.join(prefix, key))


    def __load_uge_snapshot_helper(self, dct, ti_path, content):
        if len(ti_path) == 1:
            dct[ti_path[0]] = content
        else:
            subdict = dct.setdefault(ti_path[0], dict())
            self.__load_uge_snapshot_helper(subdict, ti_path[1:], content)


