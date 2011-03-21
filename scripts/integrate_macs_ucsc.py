#!/usr/bin/env python

import getpass
import os
import sys
from optparse import OptionParser
from pypeline import Pypeline, ProcessPypeStep as PPS, PythonPypeStep as PyPS

from chipsequtil import get_org_settings

usage = "%prog <org> <stage dir> <stage url> <MACS wiggle directory>"
description = """Process a MACS wiggle directory when macs is invoked
with --wig option, convert all gzipped chromosome wiggle files to
bigWig format, copy to web staging directory <stage dir>, and create
track lines for adding to UCSC genome browser.  Requires a <org> argument
that has a path using *org_settings.py <org> ucsc_chrom_sizes* that
points to a sizes file as created by UCSC's *fetchChromSizes <org>*
tool."""

parser = OptionParser(usage=usage,description=description)
parser.add_option('--auto',dest='auto',action='store_true',help='run all steps non-interactively (for batch mode, e.g.)')

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 4 :
        parser.error('Exactly four non-option arguments required')

    organism, stage_dir, stage_url, macs_dir = args

    pipeline = Pypeline('UCSC Integration',log='ucsc_integ.log')

    steps = []

    org_settings = get_org_settings(organism)

    macs_path, macs_wiggle_path = os.path.dirname(macs_dir), os.path.basename(macs_dir)
    macs_name = macs_wiggle_path.replace('_MACS_wiggle','')
    wiggle_dir = macs_name+'_MACS_wiggle'
    bigwig_fn = macs_name+'_%s_all_chr.bw'
    d = {'wiggle_dir':macs_name+'_MACS_wiggle',
         'chrom_sizes':org_settings['ucsc_chrom_sizes'],
         'treat_bigwig_fn':macs_name+'_treat_all_chr.bw',
         'control_bigwig_fn':macs_name+'_control_all_chr.bw',
         'stage_dir':stage_dir,
         'stage_url':stage_url,
         'pwd':os.getcwd(),
        }

    # create bigWig files
    zcat_treat_call = "zcat %(wiggle_dir)s/treat/*.gz | \
                       grep -v '^track' | \
                       sed 's/\.fa//g' | \
                       wigToBigWig -clip stdin %(chrom_sizes)s \
                       %(wiggle_dir)s/treat/%(treat_bigwig_fn)s"%d
    zcat_control_call = "zcat %(wiggle_dir)s/control/*.gz |  \
                         grep -v '^track' | \
                         sed 's/\.fa//g' | \
                         wigToBigWig -clip stdin %(chrom_sizes)s \
                         %(wiggle_dir)s/control/%(control_bigwig_fn)s"%d
    steps.append(PPS('Convert wig to bigWig',[zcat_treat_call,zcat_control_call]))

    # create the staging directory
    mk_stage_dir_call = "mkdir -p %(stage_dir)s/%(wiggle_dir)s"%d
    steps.append(PPS('Create staging directory',[mk_stage_dir_call]))

    # stage bigWig files to staging directory (create links)
    stage_treat_call = "ln -fs %(pwd)s/%(wiggle_dir)s/treat/%(treat_bigwig_fn)s %(stage_dir)s/%(wiggle_dir)s/%(treat_bigwig_fn)s"%d
    stage_control_call = "ln -fs %(pwd)s/%(wiggle_dir)s/control/%(control_bigwig_fn)s %(stage_dir)s/%(wiggle_dir)s/%(control_bigwig_fn)s"%d
    steps.append(PPS('Stage bigWig files',[stage_treat_call,stage_control_call]))

    # generate track lines for treatment and control
    treat_track_d = ['track',
               'type=bigWig',
               'name="Treatment"',
               'description="%s Treatment"'%macs_name,
               'bigDataUrl=%(stage_url)s/%(wiggle_dir)s/%(treat_bigwig_fn)s'%d]
    treat_track = ' '.join(treat_track_d)

    control_track_d = ['track',
               'type=bigWig',
               'name="Control"',
               'description="%s Control"'%macs_name,
               'bigDataUrl=%(stage_url)s/%(wiggle_dir)s/%(control_bigwig_fn)s'%d]
    control_track = ' '.join(control_track_d)
    track_str = '\n'.join([treat_track,
                          control_track])

    track_fn = wiggle_dir+'_tracks.txt'
    def track_call(track_fn, track_str) :
        f = open(track_fn,'w')
        f.write(track_str+'\n')
        f.close()
    steps.append(PyPS('Generate track lines file',track_call,
                      callable_args=(track_fn,track_str))
                )

    #calls = [zcat_treat_call,
    #         zcat_control_call,
    #         mk_stage_dir_call,
    #         stage_treat_call,
    #         stage_control_call,
    #         track_call
    #         ]

    #print calls
    #steps.append(PPS('Stage Wiggle',calls))

    pipeline.add_steps(steps)
    pipeline.run(interactive=not opts.auto)
