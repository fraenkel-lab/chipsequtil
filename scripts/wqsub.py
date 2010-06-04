#!/usr/bin/env python

from __future__ import with_statement
import os
import sys
from optparse import OptionParser
from subprocess import Popen, PIPE

from chipsequtil import get_file_parts

usage = "[%prog] [options] command"
description = """Wrap the specified command into a qsub script and submit
it for execution.  Output of the command goes to file, stderr is ignored
unless specified as an option.  By default, all of the user's environment
variables are put into the script.  Environment variables that contain spaces
are not included (because SGE doesn't like them)."""
epilog = "Note: this script only works in Unix-style environments."
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('--wqsub-name',dest='wqsub_name',default='wqsub',help='job name to submit as <--wqsub-name>_<first non-whitespace chars in command> [default: %default]')
parser.add_option('--wqsub-ext',dest='wqsub_ext',default='.out',help='file extension to use for stdout files')
parser.add_option('--wqsub-keep-stderr',dest='wqsub_keep_stderr',action='store_true',help='write stderr to file (for debugging)')
parser.add_option('--wqsub-keep-script',dest='wqsub_keep_script',action='store_true',help='do not delete qsub script generated after job submission')
parser.add_option('--wqsub-no-env',dest='wqsub_no_env',action='store_true',help='do not include any local environment variables in the script')


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

    command = ' '.join(other_args)
    runscript_tmpl = """\
#!/bin/bash

#$ -N %(jobname)s
#$ -S /bin/sh
#$ -o %(stdout)s
#$ -e %(stderr)s
#$ -cwd

%(env)s

%(command)s
"""

    jobname = opts.wqsub_name+'_'+other_args[0]
    stdout_fn = jobname+opts.wqsub_ext
    stdout = os.path.abspath(stdout_fn)
    fpath,fname,fbase,fext = get_file_parts(stdout)
    stderr = '/dev/null' if not opts.wqsub_keep_stderr else os.path.abspath(os.path.join(jobname+'.err'))

    # get the user's current environment and put it into the execute script
    if opts.wqsub_no_env :
        env_str = '# local environment variables omitted'
    else :
        env_vars = ['export %s=%s'%(k,v) for k,v in os.environ.items() if v.count(' ') == 0]
        env_str = '\n'.join(env_vars)

    call_script = runscript_tmpl%{'jobname':fname,'stdout':stdout,'stderr':stderr,'command':command,'env':env_str}

    script_fn = os.path.abspath(jobname+'.script')
    with open(script_fn,'w') as f :
        f.write(call_script)
    p = Popen('qsub %s'%f.name,shell=True)
    p.wait()
    if not opts.wqsub_keep_script :
        os.remove(f.name)

