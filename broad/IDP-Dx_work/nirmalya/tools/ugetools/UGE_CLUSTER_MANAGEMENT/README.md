# UGE monitoring tools

This repository contains scripts for monitoring Broad UGE implementations.  It has been designed for UGER, but would 
probably work for any UGE implementation.  Currently there are only 2 components:

* monitor_jobs.py: A script designed to be run as a daemon that periodically checks UGE status and sends email messages
about problems.
* uge_functions.py: Tools for invoking UGE programs and parsing UGE XML output.

Folks are encouraged to suggest and/or contribute enhancements to improve the sensitivity and specificity of these tools.
Pythonistas whose style is more up to date than Alec's circa 2004 Python style are encouraged to clean up as they see fit.

**WARNING**

BITS is concerned that use of this script by many groups will cause excessive load on the qmaster.  The plan is for BITS
to run this script and address any problems reported.  If you find that BITS is not addressing problems appropriately,
you are encouraged to discuss it with them rather than run this script yourself.  It's possible that limited polling 
frequency, i.e. setting --sleep-mins and --heartbeat-mins relatively high, and monitoring the jobs of a small number
of users will be OK with BITS, but you should check with them.