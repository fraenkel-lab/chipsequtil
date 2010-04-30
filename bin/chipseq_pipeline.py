#!/usr/bin/env python

import sys, os
from pypeline import Pypeline, ProcessPypeStep as PPS
from optparse import OptionParser, OptionGroup

usage = "%prog [options] <treatment GERALD alignment filename> <control GERALD alignment filename>"
description = """1st generation ChIPSeq analysis pipeline: splits files for process distribution,
converts input files to proper format, generates alignment statistics on
treatment and control, runs MACS to find peaks, maps peaks to genes."""
epilog = """Note: it is advised to leave the --split-args, --bed-args, --qsub-args, or
--stats-args arguments unchanged unless you really know what you're doing"""

parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('--auto',dest='auto',action='store_true',help='run all steps non-interactively (for batch mode, e.g.)')
parser.add_option('--exp-name',dest='exp_name',default=os.path.basename(os.getcwd()),help='name for the experiment/pipeline, used for convenience [default: current directory name]')
parser.add_option('--split-args',dest='split_args',default='--type=count --arg=16',help='double quote wrapped arguments for split_file.py')
parser.add_option('--bed-args',dest='bed_args',default='--stdout',help='double quote wrapped arguments for gerald_to_bed.py')
parser.add_option('--stats-args',dest='stats_args',default='',help='double quote wrapped arguments for gerald_stats.py')
parser.add_option('--qsub-args',dest='qsub_args',default='',help='double quote wrapped arguments for split_qsub.py')


def get_file_parts(path) :
	path,fn = os.path.split(path)
	basename,ext = os.path.splitext(fn)
	return path,fn,basename,ext

if __name__ == '__main__' :

	# parse command line arguments
	opts, args = parser.parse_args(sys.argv[1:])

	# filenames and paths
	experiment_fn, control_fn = args[0:2]

	exp_fpath,exp_fname,exp_fbase,exp_fext = get_file_parts(experiment_fn)
	exp_wrk_dir = os.path.abspath('.exp_%s_%s'%(exp_fbase,opts.exp_name))

	cnt_fpath,cnt_fname,cnt_fbase,cnt_fext = get_file_parts(control_fn)
	cnt_wrk_dir = os.path.abspath('.cnt_%s_%s'%(cnt_fbase,opts.exp_name))

	# the pipeline
	pipeline = Pypeline()

	steps = []

	# split up files
	calls = ["mkdir %s"%exp_wrk_dir,
	         "split_file.py %s --outdir=%s %s"%(opts.split_args,exp_wrk_dir,experiment_fn),
	         "mkdir %s"%cnt_wrk_dir,
	         "split_file.py %s --outdir=%s %s"%(opts.split_args,cnt_wrk_dir,control_fn),
	        ]
	steps.append(PPS('Split files',calls))

	# convert to BED format
	exp_bed_fn = "%s_exp.bed"%exp_fbase
	cnt_bed_fn = "%s_cnt.bed"%cnt_fbase
	calls = ["split_qsub.py %s --ext=.bed gerald_to_bed.py --util-args=\"%s\" %s/*.[0-9][0-9][0-9][0-9]"%(opts.qsub_args,opts.bed_args,exp_wrk_dir),
	         "split_qsub.py %s --ext=.bed gerald_to_bed.py --util-args=\"%s\" %s/*.[0-9][0-9][0-9][0-9]"%(opts.qsub_args,opts.bed_args,cnt_wrk_dir),
	         "wait_for_qsub.py",
	         "cat %s/*.bed > %s"%(exp_wrk_dir,exp_bed_fn),
	         "cat %s/*.bed > %s"%(cnt_wrk_dir,cnt_bed_fn),
	        ]
	steps.append(PPS('Convert GERALD to BED format',calls))

	#steps.append(PPS('Helloooooooo nurse','echo Helloooooooo nurse'))
	# generate alignment statistics
	exp_stats_fn = '%s_stats.txt'%exp_fbase
	cnt_stats_fn = '%s_stats.txt'%cnt_fbase
	calls = ["split_qsub.py %s --ext=.stats gerald_stats.py --util-args=\"%s\" %s/*.[0-9][0-9][0-9][0-9]"%(opts.qsub_args,opts.stats_args,exp_wrk_dir),
	         "split_qsub.py %s --ext=.stats gerald_stats.py --util-args=\"%s\" %s/*.[0-9][0-9][0-9][0-9]"%(opts.qsub_args,opts.stats_args,cnt_wrk_dir),
	         "wait_for_qsub.py",
	         "combine_gerald_stats.py %s/*.stats > %s"%(exp_wrk_dir,exp_stats_fn),
	         "combine_gerald_stats.py %s/*.stats > %s"%(cnt_wrk_dir,cnt_stats_fn),
	        ]
	steps.append(PPS('Calculate alignment statistics',calls))

	# run macs
	calls = ["macs -t %(exp_fn)s -c %(cnt_fn)s --name=%(name)s --format=BED --tsize=35 --bw=150 --pvalue=1"%{'exp_fn':exp_bed_fn,'cnt_fn':cnt_bed_fn,'name':opts.exp_name}]
	steps.append(PPS('Run MACS',calls))

	# map peaks to genes
	# cleanup
	rm_str = "rm -f %(d)s/*.out %(d)s/*.err %(d)s/*.script %(d)s/*.stats %(d)s/*.bed"
	calls = [rm_str%{'d':exp_wrk_dir},
	         rm_str%{'d':cnt_wrk_dir}]
	steps.append(PPS('Clean up',calls))

	pipeline.add_steps(steps)
	pipeline.run(interactive=not opts.auto)
