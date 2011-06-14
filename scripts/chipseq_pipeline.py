#!/usr/bin/env python

import os
from subprocess import Popen, PIPE
import string
import sys
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP

from pypeline import Pypeline, ProcessPypeStep as PPS, PythonPypeStep as PyPS, parse_steplist
from chipsequtil import get_file_parts, get_org_settings
from chipsequtil.util import MultiLineHelpFormatter
from TAMO import MotifTools
from TAMO.MD.THEME import parser as theme_parser

usage = "%prog [options] <organism> <experiment alignment filename> [<control alignment filename>]"
description = """1st generation ChIPSeq analysis pipeline:

  - runs MACS to find peaks and sorts peaks by p-value
  - sorts peaks by pvalue and isolates top *n*
  - maps peaks to genes
  - extracts fasta files for gene peaks in experiments
  - constructs background sequences matching foreground distribution
  - runs THEME.py on input sequences w/ refinement
  - builds an infosite with stats from this analysis

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
parser.add_option('--steplist',dest='steplist',default='',help='with --auto, run specific steps')
parser.add_option('--exp-name',dest='exp_name',default=os.path.basename(os.getcwd()),help='name for the experiment/pipeline, used for convenience [default: current directory name]')
parser.add_option('--bed-args',dest='bed_args',default='--stdout --chromo-strip=.fa',help='double quote wrapped arguments for gerald_to_bed.py [default: %default]')
#parser.add_option('--stats-args',dest='stats_args',default='',help='double quote wrapped arguments for gerald_stats.py [default: %default]')
parser.add_option('--macs-exec',dest='macs_exec',default='macs14',help='the executable to use for MACS, if not an absolute path it needs to be on your shell environment path [default: %default]')
parser.add_option('--macs-args',dest='macs_args',default='--pvalue=1e-5',help='double quote wrapped arguments for macs, only changing --mfold, --tsize, --bw, and --pvalue recommended [default: %default]')
parser.add_option('--map-args',dest='map_args',default='--tss --upstream-window=10000 --downstream-window=10000',help='double quote wrapped arguments for mapping peaks to genes [default: %default]')
parser.add_option('--filter-peaks-args',dest='filter_peaks_args',default="--sort-by=pvalue --top=1000 -f 'tags>20'",help='double quote wrapped arguments for filter_macs_peaks.py [default: %default]')
parser.add_option('--filter-neg-peaks-args',dest='filter_neg_peaks_args',default="-f 'tags>20'",help='double quote wrapped arguments for filter_macs_peaks.py applied to negative peaks [default: %default]')
parser.add_option('--peaks-to-fa-args',dest='peaks_to_fa_args',default='--fixed-peak-width=200',help='double quote wrapped arguments for peaks_to_fasta.py [default: %default]')
parser.add_option('--bg-exec',dest='bg_exec',default='rejection_sample_fasta.py',help='the executable to use for generating background sequences for THEME, if not an absolute path it needs to be on your shell environment path [default: %default]')
parser.add_option('--bg-args',dest='bg_args',default='--num-seq=2.1x',help='double quote wrapped arguments for background sequence generation utility [default: %default]')
parser.add_option('--theme-args',dest='theme_args',default='--beta=0.7 --cv=5 --trials=25',help='double quote wrapped arguments for THEME.py [default: %default]')
parser.add_option('--motif-pval-cutoff',dest='motif_pval',type='float',default=1e-5,help='the p-value cutoff for sending non-refined enrichmed motifs to THEME for refinement')
parser.add_option('--parallelize',dest='parallelize',action='store_true',help='parallelize portions of the pipeline using qsub, only works from SGE execution hosts')
parser.add_option('--ucsc',dest='ucsc',action='store_true',default=False,help='perform tasks for automated integration with UCSC genome browser [default:%default]')
parser.add_option('--build-infosite-args',dest='infosite_args',default='',help='arguments to pass to build_chipseq_infosite.py [default: None]')

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
            elif opt_str == '--steplist' and not opts.auto :
                pass
            #elif opt_str in ['--stage-dir','--stage-url'] and not opts.ucsc :
            #    pass
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
        opts_strs.append('    $@')
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
    #log_fn = os.path.join(opts.exp_name+'_pipeline.log')
    pipeline = Pypeline('Analysis pipeline for %s'%opts.exp_name)

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
    macs_mfold = [x for x in opts.macs_args.split(' ') if 'mfold' in x]
    macs_mfold = macs_mfold[0].split('=',1)[1] if len(macs_mfold) >= 1 else 'DEF'

    macs_pvalue = [x for x in opts.macs_args.split(' ') if 'pvalue' in x]
    macs_pvalue = macs_pvalue[0].split('=',1)[1] if len(macs_pvalue) >= 1 else 'DEF'
    macs_name = opts.exp_name+'_mfold%s_pval%s'%(macs_mfold,macs_pvalue)

    macs_peaks_fn = macs_name+'_peaks.xls'
    macs_neg_peaks_fn = macs_name+'_negative_peaks.xls'
    macs_screen_output_fn = macs_name+'_output.txt'

    macs_d = {'exp_fn':experiment_fn,
              'cnt_flag':cnt_flag,
              'name':macs_name,
              'macs_exec':opts.macs_exec,
              'macs_args':opts.macs_args,
              'macs_out':macs_screen_output_fn,
              'gsize':org_settings['genome_size'],
              }
    calls = ["%(macs_exec)s --gsize=%(gsize)s -t %(exp_fn)s %(cnt_flag)s --name=%(name)s %(macs_args)s 2>&1 | tee %(macs_out)s"%macs_d]
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
    map_stats_fn = "%s_genes_stats.xls"%macs_name
    map_d = {'kg_ref':kg_ref,
             'kg_xref':kg_xref,
             'peaks_fn':macs_peaks_fn,
             'bed_peaks_fn':macs_name+'_peaks.bed',
             'map_fn':map_fn,
             'map_stats_fn':map_stats_fn,
             'map_args':opts.map_args
            }
    # make sure peak files don't have .fa at the end of their chromosomes
    calls = ["sed -i 's/\.fa//g' %(peaks_fn)s %(bed_peaks_fn)s"%map_d]
    c = "map_peaks_to_known_genes.py %(map_args)s --map-output=%(map_fn)s " + \
         "--detail --stats-output=%(map_stats_fn)s %(kg_ref)s %(kg_xref)s " + \
         "%(peaks_fn)s"
    calls.append(c%map_d)
    steps.append(PPS('Map peaks to genes',calls,env=os.environ))


    ############################################################################
    # filter macs peaks
    ############################################################################
    filtered_d = {'filter_peaks_args':opts.filter_peaks_args,
                  'filter_neg_peaks_args':opts.filter_neg_peaks_args,
                  'peaks_fn':macs_peaks_fn,
                  'neg_peaks_fn':macs_neg_peaks_fn
                 }
    c = "filter_macs_peaks.py --print-encoded-fn --encode-filters " \
        "%(filter_peaks_args)s %(peaks_fn)s"
    filtered_peaks_fn = Popen(c%filtered_d,shell=True,stdout=PIPE).communicate()[0]
    filtered_neg_peaks_fn = macs_name + '_negative_peak_filt.xls'
    calls = ["filter_macs_peaks.py --encode-filters %(filter_peaks_args)s %(peaks_fn)s"%filtered_d]
    if control_fn is not None :
         calls.append("filter_macs_peaks.py --encode-filters %(filter_neg_peaks_args)s %(neg_peaks_fn)s"%filtered_d)
    steps.append(PPS('Filter MACS peaks',calls,env=os.environ))


    ############################################################################
    # THEME
    ############################################################################
    # extract foreground and generate background sequences
    fg_fn = filtered_peaks_fn.replace('.xls','.fa')
    fg_d = {'opts':opts.peaks_to_fa_args,
            'organism':organism,
            'fg_fn':fg_fn,
            'peaks_fn':filtered_peaks_fn}
    calls = ["peaks_to_fasta.py %(opts)s --output=%(fg_fn)s %(organism)s %(peaks_fn)s"%fg_d]
    steps.append(PPS('Peaks to Fasta',calls,env=os.environ))

    bg_fn = "%s_bg.fa"%macs_name
    bg_d = {'opts':opts.bg_args,
            'organism':organism,
            'fg_fn':fg_fn,
            'bg_fn':bg_fn}
    calls = ["rejection_sample_fasta.py %(opts)s --output=%(bg_fn)s %(organism)s %(fg_fn)s"%bg_d]
    steps.append(PPS('Generate Background Sequences',calls,env=os.environ))

    # run THEME on fg
    theme_opts, theme_args = theme_parser.parse_args(opts.theme_args.split(' '))
    hyp_fn = org_settings['theme_hypotheses']
    markov_fn = org_settings['theme_markov']

    # run THEME w/ randomization by running each motif individuall
    # this is because TAMO.MD has a memory leak
    raw_motif_fn = '%s_motifs_beta%s_cv%s.tamo'%(macs_name,theme_opts.beta,theme_opts.cv)
    random_cv_fn = '%s_motifs_beta%s_cv%s_rand.txt'%(macs_name,theme_opts.beta,theme_opts.cv)

    # new old THEME call
    #Usage: THEME.sh [options] <FG_FASTA> <BG_FASTA> <HYP_FN> <MARKOV>
    #
    #Run old THEME version
    #
    #Options:
    #  -h, --help            show this help message and exit
    #  --hyp-indices=HYP_INDS
    #                        0-based indices of hypotheses to run [default: ALL]
    #  --no-refine           do not run with refinement
    #  --no-parallelize      do not use wqsub.py for parallelization
    #  -v, --verbose         print out the commands that are being run
    #  --dump                dump categtories to file
    #  --output-filename=OUTPUT_FN
    #                        filename to write motif results to [default:dummy.txt]
    #  --random-output=RANDOM_FN
    #                        filename to write motif results to
    #                        [default:random.txt]
    #  --motif-file=MOTIF_FN
    #                        filename to write motif results to [default:dummy.out]
    #  --beta=BETA           beta parameter to use [default:0.7]
    #  --delta=DELTA         delta parameter to use [default:0.001]
    #  --remove-common       remove common sequences from analysis
    #  --randomization       run randomization
    #  --cv=CV               number of cross validation folds [default:5]
    #  --interactive         run the script interactively

    motif_fn = '%s_motifs_beta%s_cv%s.txt'%(macs_name,theme_opts.beta,theme_opts.cv)
    theme_d = {'opts':opts.theme_args,
               'fg_fn':fg_fn,
               'bg_fn':bg_fn,
               'hyp':hyp_fn,
               'markov':markov_fn,
               'tamo_motif_fn':raw_motif_fn,
               'random_fn':random_cv_fn,
               'motif_fn':motif_fn
              }

    theme_call = "THEME.sh %(opts)s %(fg_fn)s %(bg_fn)s %(hyp)s %(markov)s " \
                 "--motif-file=%(tamo_motif_fn)s " \
                 "--random-output=%(random_fn)s " \
                 "--output-filename=%(motif_fn)s " \
                 "--randomization"

    calls = [theme_call%theme_d]
    steps.append(PPS('Run THEME',calls,env=os.environ))

    # build infosite
    calls = ['build_chipseq_infosite.py %s'%opts.infosite_args]
    steps.append(PPS('Build infosite',calls,env=os.environ))

    # cleanup
    rm_str = "rm -f %(d)s/*.out %(d)s/*.err %(d)s/*.script %(d)s/*.stats %(d)s/*.bed"
    calls = [rm_str%{'d':exp_wrk_dir}]

    if control_fn :
         calls.append(rm_str%{'d':cnt_wrk_dir})
    #steps.append(PPS('Clean up',calls,env=os.environ))

    pipeline.add_steps(steps)
    if opts.auto and opts.steplist is not None :
        steplist = parse_steplist(opts.steplist,pipeline)
    else :
        steplist = None
    pipeline.run(interactive=not opts.auto,steplist=steplist)
