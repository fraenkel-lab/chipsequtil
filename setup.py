#!/usr/bin/env python

import os
import sys
from distutils.core import setup, run_setup
from distutils.command.build_py import build_py

# convenience is king
opj = os.path.join

# make sure org_settings.cfg is in source directory
org_settings_fn = 'org_settings.cfg'
dist_settings_path = opj(os.getcwd(),'src','chipsequtil',org_settings_fn)
if not os.path.exists(dist_settings_path) :
    sys.stderr.write('WARNING: %s could not be found \
                      in distribution root directory.  org_settings.py script may \
                      not work properly.\n'%dist_settings_path)

scripts = ['scripts/chipseq_pipeline.py',
           'scripts/create_pipeline_script.py',
           'scripts/combine_gerald_stats.py',
           'scripts/gerald_stats.py',
           'scripts/gerald_to_bed.py',
           'scripts/rejection_sample_fasta.py',
           'scripts/peaks_to_fasta.py',
           'scripts/filter_macs_peaks.py',
           'scripts/sort_bed.py',
           'scripts/split_file.py',
           'scripts/split_qsub.py',
           'scripts/wqsub.py',
           'scripts/wqsub_drmaa.py',
           'scripts/wait_for_qsub.py',
           'scripts/map_peaks_to_genes.py',
           'scripts/map_peaks_to_known_genes.py',
           'scripts/join_mapped_known_genes.py',
           'scripts/compare_microarray_binding.py',
           'scripts/filter_mapped_known_genes.py',
           'scripts/filter_bed_by_position_count.py',
           'scripts/org_settings.py',
           'scripts/nibFrag.py',
           'scripts/generate_stats_doc.py',
           'scripts/probeset_to_known_gene.py',
           'scripts/compile_THEME_results.py',
           'scripts/integrate_macs_ucsc.py'
           ]

# setup and install
setup(name='chipsequtil',
      version='0.1',
      author='Adam Labadorf',
      author_email='alabadorf@gmail.com',
      package_dir={'':'src'},
      py_modules=['chipsequtil.nib','chipsequtil.util','chipsequtil.plotting','chipsequtil.sampling'],
      packages=['chipsequtil'],
      package_data={'': ['org_settings.cfg']},
      scripts=scripts,
      #cmdclass={'uninstall': uninstall},
     )
