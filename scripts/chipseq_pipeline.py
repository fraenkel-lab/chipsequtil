#!/usr/bin/env python

import os
from subprocess import Popen, PIPE
import sys
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP

from pypeline import Pypeline, ProcessPypeStep as PPS
from chipsequtil import get_file_parts, get_org_settings
from chipsequtil.util import MultiLineHelpFormatter
from TAMO2.MD.THEME import parser as theme_parser

usage = "%prog [options] <organism> <experiment BED alignment filename> [<control BED alignment filename>]"
description = """1st generation ChIPSeq analysis pipeline:

  - runs MACS to find peaks and sorts peaks by p-value
  - sorts peaks by pvalue and isolates top *n*
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
 - org_settings.py <organsim> refgene_anno_path

If running THEME:
 - org_settings.py <organism> theme_hypotheses
 - org_settings.py <organism> theme_markov

"""

epilog = """Note: it is advised to leave the --*-args arguments unchanged
unless you really know what you're doing."""

parser = OptionParser(usage=usage,description=description,epilog=epilog,formatter=MultiLineHelpFormatter())
parser.add_option('--auto',dest='auto',action='store_true',help='run all steps non-interactively (for batch mode, e.g.)')
parser.add_option('--exp-name',dest='exp_name',default=os.path.basename(os.getcwd()),help='name for the experiment/pipeline, used for convenience [default: current directory name]')
parser.add_option('--bed-args',dest='bed_args',default='--stdout --chromo-strip=.fa',help='double quote wrapped arguments for gerald_to_bed.py [default: %default]')
#parser.add_option('--stats-args',dest='stats_args',default='',help='double quote wrapped arguments for gerald_stats.py [default: %default]')
parser.add_option('--macs-exec',dest='macs_exec',default='macs14',help='the executable to use for MACS, if not an absolute path it needs to be on your shell environment path [default: %default]')
parser.add_option('--macs-args',dest='macs_args',default='--mfold=10,30 --tsize=35 --bw=150 --pvalue=1e-5',help='double quote wrapped arguments for macs, only changing --mfold, --tsize, --bw, and --pvalue recommended [default: %default]')
parser.add_option('--map-args',dest='map_args',default='--tss --upstream-window=10000 --downstream-window=10000',help='double quote wrapped arguments for mapping peaks to genes [default: %default]')
parser.add_option('--filter-peaks-args',dest='filter_peaks_args',default='--sort-by=pvalue --top=200',help='double quote wrapped arguments for filter_macs_peaks.py [default: %default]')
parser.add_option('--peaks-to-fa-args',dest='peaks_to_fa_args',default='',help='double quote wrapped arguments for peaks_to_fasta.py [default: %default]')
parser.add_option('--bg-exec',dest='bg_exec',default='rejection_sample_fasta.py',help='the executable to use for generating background sequences for THEME, if not an absolute path it needs to be on your shell environment path [default: %default]')
parser.add_option('--bg-args',dest='bg_args',default='--num-seq=2x',help='double quote wrapped arguments for background sequence generation utility [default: %default]')
parser.add_option('--theme-args',dest='theme_args',default='--beta=0.7 --cv=5',help='double quote wrapped arguments for THEME.py [default: %default]')
#parser.add_option('--parallelize',dest='parallelize',default=False,action='store_true',help='parallelize portions of the pipeline using qsub, only works from SGE execution hosts')
parser.add_option('--ucsc',dest='ucsc',action='store_true',default=False,help='perform tasks for automated integration with UCSC genome browser [default:%default]')

ucsc_group = OptionGroup(parser,"UCSC Integration Options (with --ucsc)")
ucsc_group.add_option('--stage-dir',dest='stage_dir',default='./',help='root directory where UCSC integration files should be made available [default: %default]')
ucsc_group.add_option('--stage-url',dest='stage_url',default='http://localhost/',help='URL where UCSC integration files will be made available over the web [default: %default]')
parser.add_option_group(ucsc_group)

#parallel_group = OptionGroup(parser,"Parallelization Options (with --parallelize)",description="These options are relevant to parallelization of the pipeline, functionality is in beta status until further notice")
#parallel_group.add_option('--split-args',dest='split_args',default='--type=count --arg=16',help='double quote wrapped arguments for split_file.py [default: %default]')
#parallel_group.add_option('--qsub-args',dest='qsub_args',default='--die-on-err',help='double quote wrapped arguments for split_qsub.py [default: %default]')
#parser.add_option_group(parallel_group)

parser.add_option('--print-args',dest='print_args',action='store_true',help=SUPPRESS_HELP) # secret ninja option


if __name__ == '__main__' :

    # parse command line arguments
    opts, args = parser.parse_args(sys.argv[1:])

    # stick it up here, so when we print out args it's updated
    if opts.ucsc and opts.macs_args.find('--wig') == -1 :
        opts.macs_args += " --wig"

    #  just print out all options as passed in for script generating purposes
    if opts.print_args :
        opts_strs = []
        all_opts = []
        all_opts.extend(parser.option_list)
        all_opts.extend(*[x.option_list for x in parser.option_groups])
        for opt in all_opts :
            opt_str = opt.get_opt_string()
            if opt_str in ['--help','--print-args'] :
                pass
            elif opt_str in ['--stage-dir','--stage-url'] and not opts.ucsc :
                pass
            #elif opt_str in ['--split-args','--qsub-args'] and not opts.parallelize :
            #    pass
            elif opt.action == 'store' :
                arg = str(getattr(opts,opt.dest))
                if arg.count(' ') > 0 or arg.find(' -') != -1 or arg.startswith('-') or arg.find('--') != -1 :
                    opts_strs.append('    %s="%s"'%(opt.get_opt_string(),str(getattr(opts,opt.dest))))
                else :
                    opts_strs.append('    %s=%s'%(opt.get_opt_string(),str(getattr(opts,opt.dest))))
            elif opt.action == 'store_true' and getattr(opts,opt.dest) :
                opts_strs.append('    %s'%opt.get_opt_string())
        sys.stdout.write(' \\\n'.join(opts_strs)+'\n')
        sys.exit(0)

    if len(args) < 2 :
        parser.error('Must provide two non-option arguments')

    # filenames and paths
    organism, experiment_fn = args[0:2]
    control_fn = None
    if len(args) > 2 :
        control_fn = args[2]

    org_settings = get_org_settings(organism)
    refgene_fn = org_settings['refgene_anno_path']
    kg_ref = org_settings['known_gene_anno_path']
    kg_xref = org_settings['known_gene_xref_path']

    exp_fpath,exp_fname,exp_fbase,exp_fext = get_file_parts(experiment_fn)
    exp_wrk_dir = os.path.abspath('.exp_%s_%s'%(exp_fbase,opts.exp_name))

    if control_fn :
        cnt_fpath,cnt_fname,cnt_fbase,cnt_fext = get_file_parts(control_fn)
        cnt_wrk_dir = os.path.abspath('.cnt_%s_%s'%(cnt_fbase,opts.exp_name))

    # the pipeline
    pipeline = Pypeline()

    steps = []

    #if opts.parallelize :
    #    # split up files
    #    calls = ["mkdir %s"%exp_wrk_dir,
    #             "split_file.py %s --outdir=%s %s"%(opts.split_args,exp_wrk_dir,experiment_fn),]
    #    if control_fn :
    #            calls.extend(["mkdir %s"%cnt_wrk_dir,
    #             "split_file.py %s --outdir=%s %s"%(opts.split_args,cnt_wrk_dir,control_fn),
    #            ])
    #    steps.append(PPS('Split files',calls,env=os.environ))

    ############################################################################
    # run macs
    ############################################################################
    cnt_flag = ''
    if control_fn :
        cnt_flag = '-c %s'%control_fn

    # parse macs_args so we can extract mfold and pvalue...in a rather silly way
    macs_mfold = [x for x in opts.macs_args.split(' ') if x.find('mfold') != -1]
    macs_mfold = macs_mfold[0].split('=',1)[1] if len(macs_mfold) >= 1 else 'DEF'

    macs_pvalue = [x for x in opts.macs_args.split(' ') if x.find('pvalue') != -1]
    macs_pvalue = macs_pvalue[0].split('=',1)[1] if len(macs_pvalue) >= 1 else 'DEF'
    macs_name = opts.exp_name+'_mfold%s_pval%s'%(macs_mfold,macs_pvalue)

    macs_peaks_fn = macs_name+'_peaks.xls'
    macs_screen_output_fn = macs_name+'_output.txt'

    macs_d = {'exp_fn':experiment_fn,
              'cnt_flag':cnt_flag,
              'name':macs_name,
              'macs_exec':opts.macs_exec,
              'macs_args':opts.macs_args,
              'macs_out':macs_screen_output_fn,
              'gsize':org_settings['genome_size'],
              }
    calls = ["%(macs_exec)s --gsize=%(gsize)s -t %(exp_fn)s %(cnt_flag)s --name=%(name)s --format=BED %(macs_args)s 2>&1 | tee %(macs_out)s"%macs_d]
    steps.append(PPS('Run MACS',calls,env=os.environ))


    ############################################################################
    # process and stage wiggle files
    ############################################################################
    if opts.ucsc :
        wiggle_dir = macs_name+'_MACS_wiggle'
        ucsc_d = {'org':organism,
                  'stage_dir':opts.stage_dir,
                  'stage_url':opts.stage_url,
                  'macs_dir':wiggle_dir,
                 }

        calls = ["integrate_macs_ucsc.py --auto %(org)s %(stage_dir)s %(stage_url)s %(macs_dir)s"%ucsc_d]
        steps.append(PPS("UCSC Integration",calls))


    ############################################################################
    # map peaks to genes
    ############################################################################
    map_fn = "%s_genes.txt"%macs_name
    map_stats_fn = "%s_genes_stats.txt"%macs_name
    #calls = ["map_peaks_to_genes.py %(refGene_fn)s %(peaks_fn)s --map-output=%(map_fn)s --stats-output=%(map_stats_fn)s"%{'refGene_fn':refgene_fn,'peaks_fn':macs_peaks_fn,'map_fn':map_fn,'map_stats_fn':map_stats_fn}]
    calls = ["map_peaks_to_known_genes.py %(map_args)s --map-output=%(map_fn)s --detail --stats-output=%(map_stats_fn)s %(kg_ref)s %(kg_xref)s %(peaks_fn)s"%{'kg_ref':kg_ref,'kg_xref':kg_xref,'peaks_fn':macs_peaks_fn,'map_fn':map_fn,'map_stats_fn':map_stats_fn,'map_args':opts.map_args}]
    steps.append(PPS('Map peaks to genes',calls,env=os.environ))


    ############################################################################
    # filter macs peaks
    ############################################################################
    filtered_d = {'filter_peaks_args':opts.filter_peaks_args,'peaks_fn':macs_peaks_fn}
    filtered_peaks_fn = Popen("filter_macs_peaks.py --print-encoded-fn --encode-filters %(filter_peaks_args)s %(peaks_fn)s"%filtered_d,shell=True,stdout=PIPE).communicate()[0]
    calls = ["filter_macs_peaks.py --encode-filters %(filter_peaks_args)s %(peaks_fn)s"%filtered_d]
    steps.append(PPS('Filter MACS peaks',calls,env=os.environ))


    ############################################################################
    # THEME
    ############################################################################
    # extract foreground and generate background sequences
    fg_fn = filtered_peaks_fn.replace('.xls','.fa')
    calls = ["peaks_to_fasta.py %(opts)s --output=%(fg_fn)s %(organism)s %(peaks_fn)s"%{'opts':opts.peaks_to_fa_args,'organism':organism,'fg_fn':fg_fn,'peaks_fn':filtered_peaks_fn},]
    steps.append(PPS('Peaks to Fasta',calls,env=os.environ))

    bg_fn = "%s_bg.fa"%macs_name
    calls = ["rejection_sample_fasta.py %(opts)s --output=%(bg_fn)s %(organism)s %(fg_fn)s"%{'opts':opts.bg_args,'organism':organism,'fg_fn':fg_fn,'bg_fn':bg_fn}]
    steps.append(PPS('Generate Background Sequences',calls,env=os.environ))

    # run THEME on fg
    theme_opts, theme_args = theme_parser.parse_args(opts.theme_args.split(' '))
    hyp_fn = org_settings['theme_hypotheses']
    markov_fn = org_settings['theme_markov']
    #calls = ["THEME.py %(opts)s --motif-file=%(motif_fn)s %(fg_fn)s %(bg_fn)s %(hyp)s %(markov)s"%{'opts':opts.theme_args,'motif_fn':motif_fn,'fg_fn':fg_fn,'bg_fn':bg_fn,'hyp':hyp_fn,'markov':markov_fn}]
    #steps.append(PPS('Run THEME on foreground',calls,env=os.environ))

    # run THEME randomization
    random_cv_fn = '%s_motifs_beta%s_cv%s_rand.txt'%(macs_name,theme_opts.beta,theme_opts.cv)
    raw_motif_fn = '%s_motifs_beta%s_cv%s.tamo'%(macs_name,theme_opts.beta,theme_opts.cv)
    calls = ["THEME2.py %(opts)s --motif-file=%(motif_fn)s --randomization --random-output=%(cv_fn)s %(fg_fn)s %(bg_fn)s %(hyp)s %(markov)s"%{'opts':opts.theme_args,'cv_fn':random_cv_fn,'fg_fn':fg_fn,'bg_fn':bg_fn,'hyp':hyp_fn,'markov':markov_fn,'motif_fn':raw_motif_fn}]
    steps.append(PPS('Run THEME w/ randomization',calls,env=os.environ))

    # compile THEME results
    motif_fn = '%s_motifs_beta%s_cv%s.txt'%(macs_name,theme_opts.beta,theme_opts.cv)
    calls = ["compile_THEME_results.py %(tamo_motif_fn)s %(random_fn)s > %(motif_fn)s"%{'tamo_motif_fn':raw_motif_fn,'random_fn':random_cv_fn,'motif_fn':motif_fn}]
    steps.append(PPS('Compile THEME motif results',calls,env=os.environ))

    # cleanup
    rm_str = "rm -f %(d)s/*.out %(d)s/*.err %(d)s/*.script %(d)s/*.stats %(d)s/*.bed"
    calls = [rm_str%{'d':exp_wrk_dir}]

    if control_fn :
         calls.append(rm_str%{'d':cnt_wrk_dir})
    #steps.append(PPS('Clean up',calls,env=os.environ))

    pipeline.add_steps(steps)
    pipeline.run(interactive=not opts.auto)
