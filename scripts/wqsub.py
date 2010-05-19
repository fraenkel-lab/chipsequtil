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
unless specified as an option."""
epilog = "Note: this script only works in Unix-style environments"
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('--wqsub-name',dest='wqsub_name',default='wqsub',help='job name to submit as <--wqsub-name>_<first non-whitespace chars in command> [default: %default]')
parser.add_option('--wqsub-ext',dest='wqsub_ext',default='.out',help='file extension to use for stdout files')
parser.add_option('--wqsub-keep-stderr',dest='wqsub_keep_stderr',action='store_true',help='write stderr to file (for debugging)')
parser.add_option('--wqsub-keep-scripts',dest='wqsub_keep_scripts',action='store_true',help='do not delete qsub script generated after job submission')


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) == 0 :
        parser.error('Must provide a command')

    command = ' '.join(args)
    runscript_tmpl = """
#!/bin/bash

#$ -N %(jobname)s
#$ -S /bin/sh
#$ -o %(stdout)s
#$ -e %(stderr)s
export PYTHONPATH=%(pythonpath)s:${PYTHONPATH}

%(command)s"""

    jobname = opts.wqsub_name+'_'+args[0]
    stdout_fn = jobname+opts.wqsub_ext
    stdout = os.path.abspath(stdout_fn)
    fpath,fname,fbase,fext = get_file_parts(stdout)
    stderr = '/dev/null' if not opts.wqsub_keep_stderr else os.path.join(jobname+'.err')

    call_script = runscript_tmpl%{'jobname':fname,'stdout':stdout,'stderr':stderr,'pythonpath':os.environ.get('PYTHONPATH',''),'command':command}

    script_fn = os.path.abspath(jobname+'.script')
    with open(script_fn,'w') as f :
        f.write(call_script)
    p = Popen('qsub %s'%f.name,shell=True)
    p.wait()
    if not opts.wqsub_keep_scripts :
        os.remove(f.name)

