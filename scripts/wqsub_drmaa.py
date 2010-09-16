#!/usr/bin/env python

from __future__ import with_statement
import os
import sys
from optparse import OptionParser
from subprocess import Popen, PIPE

import drmaa

from chipsequtil import get_file_parts

usage = "[%prog] [options] command"
description = """Submit *command* to a DRMAA-enabled job queueing system.
Output of the command goes to file, stderr is ignored unless specified 
as an option.  By default, all of the user's environment
variables are imported into job environment."""
epilog = "Note: this script only works in Unix-style environments."
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('--wqsub-name',dest='wqsub_name',default='wqsub',help='job name to submit as <--wqsub-name>_<first non-whitespace chars in command> [default: %default]')
parser.add_option('--wqsub-stdout',dest='wqsub_stdout',default=None,help='name of file to write stdout to (equivalent to -o argument in SGE) [default: <wqsub-name>_<command>.out]')
parser.add_option('--wqsub-stderr',dest='wqsub_stderr',default=None,help='name of file to write stderr to (equivalent to -e argument in SGE) [default: <wqsub-name>_<command>.err]')
parser.add_option('--wqsub-join',dest='wqsub_join',action='store_true',help='join stdout and stderr into file indicated by --wqsub-stdout (equivalent to -j flag in SGE)')
parser.add_option('--wqsub-no-env',dest='wqsub_no_env',action='store_true',help='do not include any local environment variables in the script')
parser.add_option('--wqsub-wait',dest='wqsub_wait',action='store_true',help='wait for job to finish executing before returning from script')


if __name__ == '__main__' :

    # get the wqsub args out first
    wqsub_args = []
    other_args = []
    for arg in sys.argv :
        if arg.count('wqsub') != 0 or arg in ['-h','--help'] :
            wqsub_args.append(arg)
        else :
            other_args.append(arg)

    opts, args = parser.parse_args(wqsub_args)

    if len(other_args) == 0 :
        parser.error('Must provide a command')

    # set up job parameters
    jobname = opts.wqsub_name+'_'+other_args[0]
    stdout_fn = jobname+'.out'
    if opts.wqsub_stdout :
        stdout_fn = opts.wqsub_stdout
    stdout = os.path.abspath(stdout_fn)

    if os.path.exists(stdout) :
        os.remove(stdout)

    stderr_fn = jobname+'.err'
    if opts.wqsub_stderr :
        stderr_fn = opts.wqsub_stderr
    stderr = os.path.abspath(stderr_fn)
    if os.path.exists(stderr) :
        os.remove(stderr)

    # drmaa job submission
    session = drmaa.Session()
    session.initialize()

    # initialize job template
    job_template = session.createJobTemplate()

    # construct DRMAA job
    command,args = other_args[0],other_args[1:]
    job_template.remoteCommand = command
    job_template.args = args
    job_template.jobName = jobname
    job_template.joinFiles = opts.wqsub_join

    # output and error paths apparently require a ':' in front
    job_template.outputPath = ':'+stdout
    job_template.errorPath = ':'+stderr

    # get the user's current environment and put it into the execute script
    if not opts.wqsub_no_env :
        job_template.jobEnvironment = os.environ

    # submit the job and wait for it
    jobid = session.runJob(job_template)

    if opts.wqsub_wait :
        # submit and wait for job to complete, keyboard interrupt aborts job
        try :

            retval = session.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)

        except KeyboardInterrupt :
            sys.stderr.write('Keyboard interrupt caught (^C), aborting')
            pass

    # clean up
    session.deleteJobTemplate(job_template)
    session.exit()
