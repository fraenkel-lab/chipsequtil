#!/usr/bin/env python

from __future__ import with_statement
import os
import re
import sys
import time
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
parser.add_option('--wqsub-drm',dest='drm',default='SGE',type='choice',choices=['SGE','TORQUE'],help='the DRM to generate scripts for [default: %default]')
parser.add_option('--wqsub-drm-arg',dest='drm_args',action='append',default=[],help='arguments to pass as parameters in the job script specific to the DRM, use multiple option flags to specify multiple parameters')
parser.add_option('--wqsub-wait',dest='wait',action='store_true',help='poll the DRM and do not return control until job is finished (only works for TORQUE)')

templates = {
'TORQUE': """\
#!/bin/bash

#PBS -N %(jobname)s
#PBS -o %(stdout)s
#PBS -e %(stderr)s
#PBS -d %(cwd)s
%(env)s
%(addnl)s

%(command)s
""",
'SGE':"""\
#!/bin/bash

#$ -N %(jobname)s
#$ -S /bin/bash
#$ -o %(stdout)s
#$ -e %(stderr)s
#$ -cwd
%(env)s
%(addnl)s

%(command)s
"""
}

drm_symb = {
'TORQUE': 'PBS',
'SGE': '$'
}

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
    runscript_tmpl = templates[opts.drm]
    # set up job parameters
    cmd_exe = os.path.basename(other_args[0])
    jobname = opts.wqsub_name+'_'+cmd_exe
    stdout_fn = jobname+opts.wqsub_ext
    stdout = os.path.abspath(stdout_fn)
    fpath,fname,fbase,fext = get_file_parts(stdout)
    stderr = os.path.abspath(os.path.join(jobname+'.err'))

    # get the user's current environment and put it into the execute script
    if opts.wqsub_no_env :
        env_str = '# local environment variables omitted'
    else :
        env_str = '#%s -V'%drm_symb[opts.drm]

    # construct the script
    addnl_params = []
    for addnl in opts.drm_args :
        addnl_params.append('#%s %s'%(drm_symb[opts.drm],addnl))
    addnl_params = '\n'.join(addnl_params)

    job_dict = {'jobname':fname,
                'stdout':stdout,
                'stderr':stderr,
                'command':command,
                'env':env_str,
                'cwd':os.getcwd(),
                'addnl':addnl_params}

    call_script = runscript_tmpl%job_dict
    # write the script to file
    script_fn = os.path.abspath(jobname+'.script')
    with open(script_fn,'w') as f :
        f.write(call_script)

    if not opts.wqsub_no_sub :
        p = Popen('qsub %s'%f.name,shell=True,stdout=PIPE)
        p.wait()
        stdout, stderr = p.communicate()
        if not opts.wqsub_keep_script :
            os.remove(f.name)
        if opts.wait :
            done = False
            print 'Waiting on job id %s'%stdout.strip()
            while not done :
                qstat_p = Popen('qstat %s'%stdout,shell=True,stdout=PIPE,stderr=PIPE)
                qstat_p.wait()
                if opts.drm == 'TORQUE' :
                    done = False if qstat_p.returncode != 153 else True
                elif opts.drm == 'SGE' :
                    done = False if qstat_p.returncode != 1 else True
                time.sleep(3) # wait three seconds because it's nice
        else :
            if opts.drm == 'TORQUE' :
                print stdout.strip()
            elif opts.drm == 'SGE' :
                qsub_output_patt = 'Your job (\d+)'
                m = re.match(qsub_output_patt,stdout.strip())
                if m is not None:
                    print m.group(1)
                    sys.exit(0)

                # might be an array job
                qsub_output_patt = 'Your job-array (\d+)\.'
                m = re.match(qsub_output_patt,stdout.strip())
                if m is not None:
                    print m.group(1)
