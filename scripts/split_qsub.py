#!/usr/bin/env python

from __future__ import with_statement
import os
import sys
from optparse import OptionParser
from subprocess import Popen, PIPE

from chipsequtil import get_file_parts

usage = "[%prog] [options] <utility> <file> [<file> <file> ...]"
description = """\
Submit a job using qsub for <utility>, each with one <file> as an argument.  Any
options specified on the command line that [%prog] cannot interpret are passed
on to the utility for each call."""
epilog = "Note: this script only works in Unix-style environments"
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('--suffix',dest='suffix',default=None,help='string to append to stdout files, e.g. <filename>_<--suffix>.<--ext> [default: <utility>]')
parser.add_option('--ext',dest='ext',default='.out',help='file extension to use for stdout files')
parser.add_option('--util-args',dest='util_args',default='',help='double quote wrapped arguments to pass to <utility>')
parser.add_option('--keep-stderr',dest='keep_stderr',action='store_true',help='capture stderr files, useful for debugging')
parser.add_option('--keep-scripts',dest='keep_scripts',action='store_true',help='do not delete qsub scripts generated after job submission')
parser.add_option('--die-on-error',dest='die_on_err',action='store_true',help='if any one of the qsub submissions returns non-zero exit status, stop executing')


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    utility, filenames = args[0], args[1:]

    # try to find the utility
    abs_utility = os.path.abspath(utility)
    if not os.path.exists(abs_utility) :
        # look on the path
        abs_utility = Popen('which %s'%utility,shell=True,stdout=PIPE,stderr=PIPE).communicate()[0].strip()
        if not os.path.exists(abs_utility) :
            raise Exception("Utility %s could not be found in the local directory or on the user's path, exiting"%utility)
            sys.exit(1)

    upath,uname,ubase,uext = get_file_parts(abs_utility)

    runscript_tmpl = """
#!/bin/bash

#$ -N %(jobname)s
#$ -S /bin/sh
#$ -o %(stdout)s
#$ -e %(stderr)s
export PYTHONPATH=%(pythonpath)s:${PYTHONPATH}

%(utility)s %(utilargs)s %(filename)s"""

    suffix = ubase if opts.suffix is None else opts.suffix
    for fn in filenames :
        abs_fn = os.path.abspath(fn)
        fpath,fname,fbase,fext = get_file_parts(abs_fn)
        stdout = os.path.join(fpath,fname+'_'+suffix+opts.ext)
        stderr = '/dev/null' if not opts.keep_stderr else os.path.join(fpath,fname+'_'+suffix+'.err')
        call_script = runscript_tmpl%{'jobname':fname,'utility':abs_utility,'filename':abs_fn,'stdout':stdout,'stderr':stderr,'utilargs':opts.util_args,'pythonpath':os.environ.get('PYTHONPATH','')}
        f = open('%s'%abs_fn+'_'+utility+'.script','w')
        f.write(call_script)
        f.close()
        p = Popen('qsub %s'%f.name,shell=True)
        p.wait()
        if not opts.keep_scripts :
            os.remove(f.name)

        if opts.die_on_err and p.returncode != 0 :
            with open(stderr,'w') as f :
                f.write('qsub returned non-zero exit code for file %s, aborting\n'%fn)
            sys.exit(1)
