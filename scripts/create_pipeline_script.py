#!/usr/bin/env python

from __future__ import with_statement
import getpass
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
utility and produces an executable run script with helpful
information on how to run it."""
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
    sys.stdout.write(white(bold(st))+'\n')


def input(st,default=None) :

    if default is None :
        default_str = ''
    else :
        default_str = ' [default: ' + default + ' ] '

    out = None
    while out is None :
        out = raw_input(white(bold(st))+default_str+white(bold(':'))+' ')
        if len(out) == 0 :
            out = default

    return out


if __name__ == '__main__' :

    TERM_ESCAPE = True

    try :

        # herro
        announce('ChIPSeq Experiment Pipeline Script Generator')
        print textwrap.fill(start_text)

        opts, args = parser.parse_args(sys.argv[1:])
        if len(args) > 0 :
            warn("Arguments were passed, but this script doesn't accept any arguments, rudely ignoring them...\n")

        ############################################################################
        # name of the experiment
        ############################################################################
        def_path = os.path.split(os.getcwd())[1]
        exp_name = input('Experiment name',def_path)
        exp_name = exp_name.replace(' ','_') # shhhhhhhh...

        ############################################################################
        # experiment and control BED file
        ############################################################################
        exp_path = input('Experiment BED alignment path')
        exp_path = exp_path.strip()
        cntrl_path = input('Control BED alignment path (leave blank for no control)','none')
        cntrl_path = cntrl_path.strip()
        if cntrl_path == 'none' :
            cntrl_path = ''

        if cntrl_path == '' :
            print 'Analysis will be run with no control'

        ############################################################################
        # organism + settings
        ############################################################################
        announce('Organism settings configuration')
        global_settings = get_global_settings()
        local_settings = get_local_settings()
        valid_org_settings = global_settings.keys() + local_settings.keys()

        org_text = """\
Below are the organism settings available on this system.  The pipeline will
use the settings for one organism (e.g. %(org)s) for the entire execution. If
you do not see a set of settings that correspond to files you need you may
add your own to %(local_org)s.  See %(glob_org)s for details.
"""

        print textwrap.fill(org_text%{'org':valid_org_settings[0],'local_org':LOCAL_SETTINGS_FN,'glob_org':GLOBAL_SETTINGS_FN},break_long_words=False)
        print

        wb('Available settings')
        # global settings
        print 'Global settings: (%s)'%GLOBAL_SETTINGS_FN
        for org, settings in global_settings.items() :
            print org
            for k,v in settings.items() :
                print ' '*4+k+": "+str(v)

        # local settings
        print 'Local settings: (%s)'%LOCAL_SETTINGS_FN
        for org, settings in local_settings.items() :
            print org
            for k,v in settings.items() :
                print ' '*4+k+": "+str(v)
        org = ''
        all_settings = {}
        all_settings.update(global_settings)
        all_settings.update(local_settings)

        while org not in valid_org_settings :
            org = input('Choose organism configuration, one of ('+','.join(valid_org_settings)+')')

            # check for the required settings
            required_settings = ['genome_dir','refgene_anno_path','theme_hypotheses','theme_markov']
            if not check_org_settings(org,required_settings) :
                warn(textwrap.fill('Selected organism settings must have the following paths defined:\n\
                     %s\n\
                     Either select another organism or define these settings in your local\
                     configuration file.'%required_settings))
                org = ''
        print

        ############################################################################
        # UCSC
        ############################################################################

        ucsc_text = """The pipeline can include a step to automatically make called
peak data available on the web for integration with UCSC genome browser."""

        print textwrap.fill(ucsc_text,break_long_words=False)

        ucsc_integrate = input('Would you like to integrate this analysis with UCSC genome browser [y/n]?','y')
        ucsc_integrate = False if ucsc_integrate == 'n' else True
        ucsc_args = ''
        if ucsc_integrate :
            stage_dir = '/nfs/antdata/web_stage/%s'%getpass.getuser()
            stage_url = 'http://fraenkel.mit.edu/stage/%s'%getpass.getuser()
            ucsc_args = ['--ucsc',
                         '--stage-dir=%s'%stage_dir,
                         '--stage-url=%s'%stage_url]
            ucsc_args = ' '.join(ucsc_args)

        # TODO - consider letting user set these on script creation time
        # any utility specific arguments?
        #  - MACS
        #  - THEME


        ############################################################################
        # various pipeline parameters
        ############################################################################
        pipeline_args = {}

        # --macs-args
        macs_args = ['--mfold=10,30','--tsize=35','--bw=150']
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
        pipeline_args['--filter-peaks-args'] = ' '.join(filt_args)

        top = ''
        while not re.search('^\d+$',top) and top != 'ALL' :
            top = input('How many peak sequences should be used for motif discovery when sorted by p-value [<# peaks>/ALL]','200')
        if top != 'ALL' :
            filt_args.append('--top=%s'%top)
        pipeline_args['--filter-peaks-args'] = ' '.join(filt_args)


        ############################################################################
        # done with input, creating script and other stuff
        ############################################################################
        # if the experiment and control files are in a different directory,
        # create symlinks for them
        exp_dir,exp_fn = os.path.split(os.path.abspath(exp_path))
        if exp_dir != os.getcwd() :
            wb('Creating symlink for experiment file...')
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
                wb('Creating symlink for control file...')
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
        def_args = Popen('chipseq_pipeline.py --exp-name=%s %s %s --print-args'%(exp_name,ucsc_args,pipeline_args),shell=True,stdout=PIPE,stderr=PIPE).communicate()[0]

        wb('Creating script...')
        script_fn = '%s_pipeline.sh'%exp_name
        with open(script_fn,'w') as script_f :
            script_f.write(script_template%{'exp_path':exp_fn,'cnt_path':cntrl_fn,'organism':org,'exp_name':exp_name,'def_args':def_args})
            os.chmod(script_f.name,stat.S_IRWXU|stat.S_IRWXG|stat.S_IROTH)

        print end_text%{'script_fn':script_fn}

    except KeyboardInterrupt :
        sys.stderr.write('\n')
        error('Script creation interrupted, aborting')
