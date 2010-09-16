#!/usr/bin/env python

import os
import sys
from optparse import OptionParser, OptionGroup

from pypeline import Pypeline, ProcessPypeStep as PPS
from chipsequtil import get_file_parts, get_org_settings
from chipsequtil.util import MultiLineHelpFormatter

usage = "%prog [options] <organism> <experiment GERALD alignment filename> [<control GERALD alignment filename>]"
description = """1st generation ChIPSeq analysis pipeline:

  - converts Illumina GERALD alignment files to BED format
  - calculates statistics on input alignments
  - runs MACS to find peaks
  - maps peaks to genes
  - extracts fasta files for gene peaks in experiments
  - constructs background sequences matching foreground distribution
  - runs THEME.py on input sequences
  - runs THEME.py randomization
  - creates documentation on entire pipeline run

Control input file is optional.  *organism* argument is passed to the
*org_settings.py* command to specify organism specific parameters, ensure
that the following commands return valid paths:

If running MACS:
 - org_settings.py <organism> genome_size
 - org_settings.py <organism> genome_dir
 - org_settings.py <organsim> annotation_path

If running THEME:
 - org_settings.py <organism> theme_hypotheses
 - org_settings.py <organism> theme_markov

"""

epilog = """Note: it is advised to leave the --*-args arguments unchanged
unless you really know what you're doing."""

parser = OptionParser(usage=usage,description=description,epilog=epilog,formatter=MultiLineHelpFormatter())
parser.add_option('--auto',dest='auto',action='store_true',help='run all steps non-interactively (for batch mode, e.g.)')
parser.add_option('--exp-name',dest='exp_name',default=os.path.basename(os.getcwd()),help='name for the experiment/pipeline, used for convenience [default: current directory name]')
parser.add_option('--split-args',dest='split_args',default='--type=count --arg=16',help='double quote wrapped arguments for split_file.py [default: %default]')
parser.add_option('--bed-args',dest='bed_args',default='--stdout --chromo-strip=.fa',help='double quote wrapped arguments for gerald_to_bed.py [default: %default]')
parser.add_option('--stats-args',dest='stats_args',default='',help='double quote wrapped arguments for gerald_stats.py [default: %default]')
parser.add_option('--qsub-args',dest='qsub_args',default='--die-on-err',help='double quote wrapped arguments for split_qsub.py [default: %default]')
parser.add_option('--macs-args',dest='macs_args',default='--mfold=10 --tsize=35 --bw=150 --pvalue=1e-5',help='double quote wrapped arguments for macs, only changing --mfold, --tsize, --bw, and --pvalue recommended [default: %default]')
parser.add_option('--pk-to-fa-args',dest='pk_to_fa_args',default='--bg-type=rej_samp',help='double quote wrapped arguments for peaks_to_fasta.py [default: %default]')
parser.add_option('--theme-args',dest='theme_args',default='--beta=0.7 --cv=5',help='double quote wrapped arguments for THEME.py [default: %default]')


if __name__ == '__main__' :

    # parse command line arguments
    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 3 :
        parser.error('Must provide two non-option arguments')

    # filenames and paths
    organism, experiment_fn, control_fn = args[0:3]
    control_fn = None
    if len(args) > 3 :
        control_fn = args[2]

    org_settings = get_org_settings(organism)
    refseq_fn = org_settings['annotation_path']

    exp_fpath,exp_fname,exp_fbase,exp_fext = get_file_parts(experiment_fn)
    exp_wrk_dir = os.path.abspath('.exp_%s_%s'%(exp_fbase,opts.exp_name))

    if control_fn :
        cnt_fpath,cnt_fname,cnt_fbase,cnt_fext = get_file_parts(control_fn)
        cnt_wrk_dir = os.path.abspath('.cnt_%s_%s'%(cnt_fbase,opts.exp_name))

    # the pipeline
    pipeline = Pypeline()

    steps = []

    # split up files
    calls = ["mkdir %s"%exp_wrk_dir,
             "split_file.py %s --outdir=%s %s"%(opts.split_args,exp_wrk_dir,experiment_fn),]
    if control_fn :
            calls.extend(["mkdir %s"%cnt_wrk_dir,
             "split_file.py %s --outdir=%s %s"%(opts.split_args,cnt_wrk_dir,control_fn),
            ])
    steps.append(PPS('Split files',calls,env=os.environ))

    # convert to BED format
    exp_bed_fn = "%s_exp.bed"%exp_fbase
    calls = ["split_qsub.py %s --ext=.bed gerald_to_bed.py --util-args=\"%s\" %s/*.[0-9][0-9][0-9][0-9]"%(opts.qsub_args,opts.bed_args,exp_wrk_dir),
             "wait_for_qsub.py",
             "cat %s/*.bed > %s"%(exp_wrk_dir,exp_bed_fn),
            ]

    if control_fn :
        cnt_bed_fn = "%s_cnt.bed"%cnt_fbase
        calls.extend(["split_qsub.py %s --ext=.bed gerald_to_bed.py --util-args=\"%s\" %s/*.[0-9][0-9][0-9][0-9]"%(opts.qsub_args,opts.bed_args,cnt_wrk_dir),
                      "wait_for_qsub.py",
                      "cat %s/*.bed > %s"%(cnt_wrk_dir,cnt_bed_fn),
                     ])

    steps.append(PPS('Convert GERALD to BED format',calls,env=os.environ))

    #steps.append(PPS('Helloooooooo nurse','echo Helloooooooo nurse'))
    # generate alignment statistics
    exp_stats_fn = '%s_stats.txt'%exp_fbase
    calls = ["split_qsub.py %s --ext=.stats gerald_stats.py --util-args=\"%s\" %s/*.[0-9][0-9][0-9][0-9]"%(opts.qsub_args,opts.stats_args,exp_wrk_dir),
             "wait_for_qsub.py",
             "combine_gerald_stats.py %s/*.stats > %s"%(exp_wrk_dir,exp_stats_fn),
            ]

    if control_fn :
        cnt_stats_fn = '%s_stats.txt'%cnt_fbase
        calls.extend(["split_qsub.py %s --ext=.stats gerald_stats.py --util-args=\"%s\" %s/*.[0-9][0-9][0-9][0-9]"%(opts.qsub_args,opts.stats_args,cnt_wrk_dir),
                 "wait_for_qsub.py",
                 "combine_gerald_stats.py %s/*.stats > %s"%(cnt_wrk_dir,cnt_stats_fn),
                ])
    steps.append(PPS('Calculate alignment statistics',calls,env=os.environ))

    # run macs
    cnt_flag = ''
    if control_fn :
        cnt_flag = '-c %s'cnt_bed_fn

    macs_d = {'exp_fn':exp_bed_fn,
              'cnt_flag':cnt_flag,
              'name':opts.exp_name,
              'macs_args':opts.macs_args,
              'gsize':org_settings['genome_size'],
              }
    calls = ["macs --gsize=%(gsize)s -t %(exp_fn)s %(cnt_flag)s --name=%(name)s --format=BED %(macs_args)s"%macs_d]
    steps.append(PPS('Run MACS',calls,env=os.environ))

    # map peaks to genes
    peaks_fn = "%s_peaks.bed"%opts.exp_name
    map_fn = "%s_genes.txt"%opts.exp_name
    map_stats_fn = "%s_genes_stats.txt"%opts.exp_name
    calls = ["map_peaks_to_genes.py --peaks-format=BED %(refGene_fn)s %(peaks_fn)s --map-output=%(map_fn)s --stats-output=%(map_stats_fn)s"%{'refGene_fn':refseq_fn,'peaks_fn':peaks_fn,'map_fn':map_fn,'map_stats_fn':map_stats_fn}]
    steps.append(PPS('Map peaks to genes',calls,env=os.environ))

    # THEME
    # extract foreground and generate background sequences
    fg_fn = "%s_peaks.fa"%opts.exp_name
    bg_fn = "%s_bg.fa"%opts.exp_name
    nib_dir = org_settings['genome_dir']
    calls = ["peaks_to_fasta.py %(opts)s --output=%(fg_fn)s --bg-fn=%(bg_fn)s %(organism)s %(peaks_fn)s"%{'opts':opts.pk_to_fa_args,'organism':organism,'fg_fn':fg_fn,'bg_fn':bg_fn,'peaks_fn':peaks_fn}]
    steps.append(PPS('Peaks to Fasta',calls,env=os.environ))

    # run THEME on fg
    motif_fn = '%s_motifs.txt'%opts.exp_name
    hyp_fn = org_settings['theme_hypotheses']
    markov_fn = org_settings['theme_markov']
    calls = ["THEME.py %(opts)s --motif-file=%(motif_fn)s %(fg_fn)s %(bg_fn)s %(hyp)s %(markov)s"%{'opts':opts.theme_args,'motif_fn':motif_fn,'fg_fn':fg_fn,'bg_fn':bg_fn,'hyp':hyp_fn,'markov':markov_fn}]
    steps.append(PPS('Run THEME on foreground',calls,env=os.environ))

    # run THEME randomization
    random_motif_fn = '%s_motifs_rand.txt'%opts.exp_name
    calls = ["THEME.py %(opts)s --randomization --motif-file=%(motif_fn)s %(fg_fn)s %(bg_fn)s %(hyp)s %(markov)s"%{'opts':opts.theme_args,'motif_fn':random_motif_fn,'fg_fn':fg_fn,'bg_fn':bg_fn,'hyp':hyp_fn,'markov':markov_fn}]
    steps.append(PPS('Run THEME randomization',calls,env=os.environ))

    # cleanup
    rm_str = "rm -f %(d)s/*.out %(d)s/*.err %(d)s/*.script %(d)s/*.stats %(d)s/*.bed"
    calls = [rm_str%{'d':exp_wrk_dir},
             rm_str%{'d':cnt_wrk_dir}]
    steps.append(PPS('Clean up',calls,env=os.environ))

    pipeline.add_steps(steps)
    pipeline.run(interactive=not opts.auto)
