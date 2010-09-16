#!/usr/bin/env python

from __future__ import with_statement
import os
import sys
from optparse import OptionParser
from subprocess import Popen, PIPE

from chipsequtil import get_file_parts

usage = "[%prog] [options] command"
description = """Wrap the specified command into a qsub script and submit it
for execution. Script captures both stdout and stderr to the current directory.
By default, all of the user's environment variables are put into the script
(compatible with SGE only ATM)."""
epilog = "Note: this script only works in Unix-style environments."
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('--wqsub-name',dest='wqsub_name',default='wqsub',help='job name to submit as <--wqsub-name>_<first non-whitespace chars in command> [default: %default]')
parser.add_option('--wqsub-ext',dest='wqsub_ext',default='.out',help='file extension to use for stdout files')
parser.add_option('--wqsub-keep-script',dest='wqsub_keep_script',action='store_true',help='do not delete qsub script generated after job submission')
parser.add_option('--wqsub-no-env',dest='wqsub_no_env',action='store_true',help='do not include any local environment variables in the script')
parser.add_option('--wqsub-no-submit',dest='wqsub_no_sub',action='store_true',help='create script but do not submit job (useful for generating scripts)')


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
#$ -S /bin/bash
#$ -o %(stdout)s
#$ -e %(stderr)s
#$ -cwd
%(env)s

%(command)s
"""

    # set up job parameters
    jobname = opts.wqsub_name+'_'+other_args[0]
    stdout_fn = jobname+opts.wqsub_ext
    stdout = os.path.abspath(stdout_fn)
    fpath,fname,fbase,fext = get_file_parts(stdout)
    stderr = os.path.abspath(os.path.join(jobname+'.err'))

    # get the user's current environment and put it into the execute script
    if opts.wqsub_no_env :
        env_str = '# local environment variables omitted'
    else :
        #env_vars = ['export %s=%s'%(k,v) for k,v in os.environ.items() if v.count(' ') == 0]
        #env_str = '\n'.join(env_vars)
        env_str = '#$ -V'

    # construct the script
    call_script = runscript_tmpl%{'jobname':fname,'stdout':stdout,'stderr':stderr,'command':command,'env':env_str}

    # write the script to file
    script_fn = os.path.abspath(jobname+'.script')
    with open(script_fn,'w') as f :
        f.write(call_script)

    if not opts.wqsub_no_sub :
        p = Popen('qsub %s'%f.name,shell=True)
        p.wait()
        if not opts.wqsub_keep_script :
            os.remove(f.name)

