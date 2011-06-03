#!/usr/bin/env python

from __future__ import with_statement
import getpass
import json
import os
import textwrap

try:
    import readline
    import glob
    readline.parse_and_bind("tab: complete")
    readline.set_completer_delims('')

    comp_states = {}
    def basic_complete_file(text,state) :
        #if text.strip() == '' :
        #    text = './'
        options = dict([(i,p) for i,p in enumerate(glob.glob(text+'*'))])
        return options.get(state,None)

    readline.set_completer(basic_complete_file)

except ImportError:
    print "Module readline not available."

import re
import stat
import sys
from optparse import OptionParser
from subprocess import Popen, PIPE

import chipsequtil
from chipsequtil import get_global_settings, get_local_settings, check_org_settings, GLOBAL_SETTINGS_FN, LOCAL_SETTINGS_FN
from terminalcontroller import TERM_ESCAPE, announce, warn, error, white, bold

usage = "%prog"
description = """Script for creating a custom run script for
ChIPSeq/DNAse hypersensitivity experiments.  User is asked for
paths and settings required for ChIPSeq analysis using the *chipseq_pipeline.py*
utility and produces an executable run script with helpful information on how to
run it.  Also creates a JSON formatted file containing all the parameters for
this pipeline run."""
epilog = "Note: this script only works in Unix-style environments"
parser = OptionParser(usage=usage,description=description,epilog=epilog)


script_template = """\
#!/bin/bash

# required parameters for the pipeline
ORG=%(organism)s
EXP_FN=%(exp_path)s
CNT_FN=%(cnt_path)s

# chipseq_pipeline.py is the main workhorse of this analysis
# you may change any of the arguments below from their defaults

chipseq_pipeline.py $ORG $EXP_FN $CNT_FN \\
%(def_args)s
"""

start_text = """\
This is an interactive script that creates an executable script to use for
ChIPSeq analyses. When prompted for experiment and control files, tab
completion is available a la bash or tcsh shells. Press Ctrl-C at any time to
quit.
"""

end_text = """The script %(script_fn)s has been created to run this pipeline. \
The script can now be run with:

$> ./%(script_fn)s

Have a nice day."""



def wb(st) :
    sys.stdout.write(white(bold(st)))


def input(st,default=None) :

    if default is None :
        default_str = ''
    else :
        default_str = ' [default: ' + default + ' ] '

    out = None
    while out is None :
        out = raw_input(white(bold(st))+default_str+white(bold(':'))+' \n')
        if len(out) == 0 :
            out = default

    return out


if __name__ == '__main__' :

    TERM_ESCAPE = True

    try :

        pipeline_args = {}

        # herro
        announce('ChIPSeq Experiment Pipeline Script Generator')
        print textwrap.fill(start_text)

        opts, args = parser.parse_args(sys.argv[1:])
        if len(args) > 0 :
            warn("Arguments were passed, but this script doesn't accept any arguments, rudely ignoring them...\n")

        # this dictionary will be used to generate a JSON formatted file with
        # all the relevant settings for the pipeline
        json_dict = {}

        ############################################################################
        # name of the experiment
        ############################################################################
        def_path = os.path.basename(os.getcwd())
        exp_name = input('Experiment name',def_path)
        exp_name = exp_name.replace(' ','_') # shhhhhhhh...

        json_dict['experiment name'] = exp_name
        json_dict['analysis path'] = os.getcwd()

        ############################################################################
        # experiment and control file
        ############################################################################
        align_text = "The pipeline can accept either BED, BOWTIE, SAM, or " \
        "ELANDEXPORT formatted alignment files. SAM is the default " \
        "format of files provided by the BMC pipeline.  Both experiment " \
        "and control files must have the same format."
        print textwrap.fill(align_text)

        align_fmt = input("Which format are the alignment files in?",'SAM')
        exp_path = input('Experiment alignment path')
        exp_path = exp_path.strip()

        lims_exp_url = input('Experiment LIMS sample URL, if applicable','none')
        lims_exp_url = lims_exp_url.strip()

        cntrl_path = input('Control alignment path (leave blank for no control)','none')
        cntrl_path = cntrl_path.strip()

        lims_cntrl_url = input('Control LIMS sample URL, if applicable','none')
        lims_cntrl_url = lims_cntrl_url.strip()

        if cntrl_path == 'none' :
            cntrl_path = ''

        if cntrl_path == '' :
            print 'Analysis will be run with no control'

        json_dict['experiment path'] = os.path.realpath(exp_path)
        json_dict['experiment lims url'] = lims_exp_url
        json_dict['control path'] = os.path.realpath(cntrl_path) if cntrl_path != '' else 'none'
        json_dict['control lims url'] = lims_cntrl_url

        ############################################################################
        # organism + settings
        ############################################################################
        announce('Organism settings configuration')
        global_settings = get_global_settings()
        local_settings = get_local_settings()
        valid_org_settings = global_settings.keys() + local_settings.keys()
        valid_org_settings.sort()

        org_text = """\
Below are the organism settings available on this system.  The pipeline will
use the settings for one organism (e.g. %(org)s) for the entire execution. If
you do not see a set of settings that correspond to files you need you may
add your own to %(local_org)s.  See %(glob_org)s for details.
"""

        print textwrap.fill(org_text%{'org':valid_org_settings[0],'local_org':LOCAL_SETTINGS_FN,'glob_org':GLOBAL_SETTINGS_FN},break_long_words=False)
        print

        wb('Available settings\n')
        # global settings
        print 'Global settings: (%s)'%GLOBAL_SETTINGS_FN
        org_sets = [(k,global_settings[k]) for k in sorted(global_settings.keys())]
        for org, settings in org_sets :
            wb(org.ljust(8))
            print ':', settings.get('description','No description')
            #for k,v in settings.items() :
            #    print ' '*4+k+": "+str(v)

        # local settings
        print 'Local settings: (%s)'%LOCAL_SETTINGS_FN
        org_sets = [(k,local_settings[k]) for k in sorted(local_settings.keys())]
        for org, settings in org_sets :
            wb(org.ljust(8))
            print ':', settings.get('description','No description')
            #for k,v in settings.items() :
            #    print ' '*4+k+": "+str(v)
        org = ''
        all_settings = {}
        all_settings.update(global_settings)
        all_settings.update(local_settings)

        while org not in valid_org_settings :
            org = input('Choose organism configuration, one of ('+','.join(valid_org_settings)+')')

            # check for the required settings
            required_settings = ['description','genome_dir','refgene_anno_path','theme_hypotheses','theme_markov']
            if not check_org_settings(org,required_settings) :
                warn(textwrap.fill('Selected organism settings must have the following settings defined:\n\
                     %s\n\
                     Either select another organism or define these settings in your local\
                     configuration file.'%required_settings))
                org = ''
        print

        json_dict['org'] = org

        ############################################################################
        # UCSC
        ############################################################################

        ucsc_text = """The pipeline can include a step to automatically make called
peak data available on the web for integration with UCSC genome browser."""

        print textwrap.fill(ucsc_text,break_long_words=False)

        ucsc_integrate = input('Would you like to integrate this analysis with UCSC genome browser [y/n]?','y')
        ucsc_integrate = False if ucsc_integrate == 'n' else True
        ucsc_args = ''
        stage_dir = '/nfs/antdata/web_stage/%s'%getpass.getuser()
        stage_url = 'http://fraenkel.mit.edu/stage/%s'%getpass.getuser()
        if ucsc_integrate :
            ucsc_args = ['--ucsc']
            ucsc_args = ' '.join(ucsc_args)

        pipeline_args['--stage-dir'] = stage_dir
        pipeline_args['--stage-url'] = stage_url

        json_dict['stage dir'] = stage_dir
        json_dict['stage url'] = stage_url

        # TODO - consider letting user set these on script creation time
        # any utility specific arguments?
        #  - MACS
        #  - THEME


        ############################################################################
        # various pipeline parameters
        ############################################################################

        # --macs-args
        macs_args = ['--mfold=10,30','--format=%s'%align_fmt]
        pval = ''
        while not re.search('^\de-\d+$',pval) :
            pval = input('What p-value should MACS use as a cutoff?','1e-5')
        macs_args.append('--pvalue=%s'%pval)
        pipeline_args['--macs-args'] = ' '.join(macs_args)

        # --map-args
        map_args = []
        tss = ''
        while tss.upper() not in ('TSS','GENE') :
            tss = input('Should gene mapping be made in relation to transcription start site or full gene coordinates [TSS/gene]?','TSS')
        if tss == 'TSS' :
            map_args.append('--tss')

        window = ''
        while not re.search('^\d+,\d+$',window) :
            window = input('What window would you like to use for mapping peaks to genes (upstream bases,downstream bases)?','10000,10000')
        upstr, downstr = window.split(',')
        map_args.extend(['--upstream-window=%s'%upstr,'--downstream-window=%s'%downstr])
        pipeline_args['--map-args'] = ' '.join(map_args)

        # --filter-peaks-args
        filt_args =  ['--sort-by=pvalue']
        fdr = ''
        while not re.search('^\d+(\.\d+)?',fdr) and fdr != 'none' :
            fdr = input('What FDR cutoff should be used, in %?','none')
        if fdr != 'none' :
            filt_args.append("--filter='fdr<%s'"%fdr)

        top = ''
        while not re.search('^\d+$',top) and top != 'ALL' :
            top = input('How many peak sequences should be used for motif discovery when sorted by p-value [<# peaks>/ALL]','1000')
        if top != 'ALL' :
            filt_args.append('--top=%s'%top)

        # tag filter for both pos and neg peaks
        tags = ''
        filt_neg_args = []
        while not re.search('^\d+$',tags) and tags != 'ALL' :
            tags = input('What tag count cutoff should be used as a minimum for positive and negative peaks? [<# peaks>/None]','20')
        if tags != 'None' :
            filt_args.append("--filter='tags>%s'"%tags)
            filt_neg_args.append("--filter='tags>%s'"%tags)
        pipeline_args['--filter-peaks-args'] = ' '.join(filt_args)
        pipeline_args['--filter-neg-peaks-args'] = ' '.join(filt_neg_args)

        # --peaks-to-fa-args
        peaks_to_fa_args = []
        width = ''
        while not re.search('^\d+$',width) and width != 'NA' :
            width = input('What width around peak summit should be used for motif analysis (NA to use entire peak)? [<# bases>/NA]','200')
        if width != 'NA' :
            peaks_to_fa_args.append('--fixed-peak-width=%s'%width)
        else :
            width = 'none'
        pipeline_args['--peaks-to-fa-args'] = ' '.join(peaks_to_fa_args)

        # --parallelize
        parallel = input('Use cluster parallelization [y/n]?','y')
        parallel = '--parallelize' if parallel.lower() != 'n' else ''

        # each user-specified argument gets its own key
        json_dict['format'] = align_fmt
        json_dict['mapping type'] = tss
        json_dict['mapping window'] = (upstr,downstr)
        json_dict['FDR filter'] = fdr
        json_dict['peaks used by THEME'] = top
        json_dict['fixed peak width'] = width
        json_dict['parallelize'] = parallel != ''
        json_dict['peak tag count filter'] = tags

        # put all the command line utility args in json_dict as its own dict
        json_dict['pipeline args'] = pipeline_args

        ############################################################################
        # done with input, creating script and other stuff
        ############################################################################
        # if the experiment and control files are in a different directory,
        # create symlinks for them
        exp_dir,exp_fn = os.path.split(os.path.abspath(exp_path))
        if exp_dir != os.getcwd() :
            wb('Creating symlink for experiment file...\n')
            if os.path.exists(exp_fn) :
                if os.path.realpath(exp_fn) != os.path.abspath(exp_path) : # existing symlink  doesn't point to the same file, prompt to overwrite
                    ans = raw_input('Symlink %s in current directory points to %s but you asked for %s, overwrite symbolic link? y/[n] '%(exp_fn,os.path.realpath(exp_fn),os.path.abspath(exp_path)))
                    if ans == 'y' :
                        os.remove(exp_fn)
                        exp_fn = 'exp_'+exp_fn
                        os.symlink(exp_path,exp_fn)
            else :
                exp_fn = 'exp_'+exp_fn
                os.symlink(exp_path,exp_fn)

        if cntrl_path != '' :
            cntrl_dir,cntrl_fn = os.path.split(os.path.abspath(cntrl_path))
            if cntrl_dir != os.getcwd() :
                wb('Creating symlink for control file...\n')
                if os.path.exists(cntrl_fn) :
                    if os.path.realpath(cntrl_fn) != os.path.abspath(cntrl_path) : # existing symlink  doesn't point to the same file, prompt to overwrite
                        ans = raw_input('Symlink %s in current directory points to %s but you asked for %s, overwrite symbolic link? y/[n] '%(cntrl_fn,os.path.realpath(cntrl_fn),os.path.abspath(cntrl_path)))
                        if ans == 'y' :
                            os.remove(cntrl_fn)
                            cntrl_fn = 'cntrl_'+cntrl_fn
                            os.symlink(cntrl_path,cntrl_fn)
                else :
                    cntrl_fn = 'cntrl_'+cntrl_fn
                    os.symlink(cntrl_path,cntrl_fn)
        else :
            cntrl_fn = ''

        # get default chipseq_pipeline.py args
        pipeline_args = ' '.join(['%s="%s"'%(k,v) for k,v in pipeline_args.items()])
        print 'chipseq_pipeline.py --exp-name=%s %s %s --print-args'%(exp_name,ucsc_args,pipeline_args)
        def_args = Popen('chipseq_pipeline.py --exp-name=%s %s %s %s --print-args'%(exp_name,ucsc_args,parallel,pipeline_args),shell=True,stdout=PIPE,stderr=PIPE).communicate()[0]

        wb('Creating script...\n')
        script_fn = '%s_pipeline.sh'%exp_name
        with open(script_fn,'w') as script_f :
            script_f.write(script_template%{'exp_path':exp_fn,'cnt_path':cntrl_fn,'organism':org,'exp_name':exp_name,'def_args':def_args})
            os.chmod(script_f.name,stat.S_IRWXU|stat.S_IRWXG|stat.S_IROTH)

        print end_text%{'script_fn':script_fn}

        wb('Creating parameter file...\n')
        json_fn = '%s_params.json'%exp_name
        with open(json_fn,'w') as json_f :
            json.dump(json_dict,json_f,indent=4)

    except KeyboardInterrupt :
        sys.stderr.write('\n')
        error('Script creation interrupted, aborting')
