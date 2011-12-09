#!/usr/bin/env python

import glob
import signal
import time
from subprocess import Popen, PIPE
from textwrap import TextWrapper

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm

signal.signal(signal.SIGALRM, alarm_handler)

scripts = [#'../scripts/build_chipseq_infosite.py',
           '../scripts/chipseq_pipeline.py',
           #'../scripts/combine_gerald_stats.py',
           #'../scripts/compare_microarray_binding.py',
           '../scripts/create_pipeline_script.py',
           '../scripts/extract_promoters.py',
           '../scripts/filter_bed_by_position_count.py',
           '../scripts/filter_macs_peaks.py',
           '../scripts/filter_gps_peaks.py',
           '../scripts/filter_mapped_known_genes.py',
           #'../scripts/generate_stats_doc.py',
           '../scripts/gerald_stats.py',
           '../scripts/gerald_to_bed.py',
           #'../scripts/integrate_macs_ucsc.py',
           '../scripts/join_mapped_known_genes.py',
           '../scripts/map_intervals.py',
           '../scripts/map_peaks_to_genes.py',
           '../scripts/map_peaks_to_known_genes.py',
           '../scripts/motif_scan.py',
           '../scripts/nibFrag.py',
           '../scripts/org_settings.py',
           '../scripts/peaks_to_fasta.py',
           '../scripts/plot_pos_vs_neg_peaks.py',
           '../scripts/plot_peak_loc_dist.py',
           #'../scripts/probeset_to_known_gene.py',
           '../scripts/rejection_sample_fasta.py',
           '../scripts/sort_bed.py',
           #'../scripts/split_file.py',
           #'../scripts/split_qsub.py',
           #'../scripts/THEME.sh',
           #'../scripts/wait_for_qsub.py',
           '../scripts/wait_for_jobid.py',
           '../scripts/wqsub.py',
           '../scripts/wqsub_drmaa.py',
           ]

if __name__ == '__main__' :

    tw = TextWrapper(initial_indent="   ",subsequent_indent="   ")
    script_help_out = ''
    refs = ''
    for script in scripts :
        cmd = 'python %s -h'%script
        p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE)

        stdout, stderr = None, None
        signal.alarm(3)  # 3 seconds
        try:
            stdout, stderr = p.communicate()
            signal.alarm(0)  # reset the alarm
        except Alarm:
            pass

        script_str = script.replace('../scripts/','')


        refs += '  - :ref:`%(script_str)s <%(script_str)s>`\n'%{'script_str':script_str}
        script_help_out += '.. _%s:\n\n'%script_str
        script_help_out += '%s::\n\n'%script_str
        if stderr is None :
            script_help_out += tw.fill('empty docstring\n')
        else :
            script_help_out += '\n'.join(['   '+x for x in stdout.split('\n')])
            script_help_out += '\n'.join(['   '+x for x in stderr.split('\n')])
        script_help_out += '\n\n'
        script_help_out += ':ref:`top <top>`\n\n'

    rst_str = """\
Illumina pipeline script reference
==================================

The following is the output of the scripts provided by this package when invoked
on the command line with *-h*.

.. _top:

Scripts:
%(refs)s

%(script_help_out)s
"""%{'refs':refs,'script_help_out':script_help_out}

    print rst_str
