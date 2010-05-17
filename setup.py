#!/usr/bin/env python
from distutils.core import setup

setup(name='chipseq',
      version='0.1',
      py_modules=['chipsequtil'],
      scripts= ['scripts/chipseq_pipeline.py',
                'scripts/combine_gerald_stats.py',
                'scripts/gerald_stats.py',
                'scripts/gerald_to_bed.py',
                'scripts/sort_bed.py',
                'scripts/split_file.py',
                'scripts/split_qsub.py',
                'scripts/wait_for_qsub.py',
                'scripts/map_peaks_to_genes.py',
                'scripts/nibFrag.py']
      package_dir={'chipsequtil':'src/chipsequtil'}
      package_data={'chipsequtil': ['org_settings.cfg']}
     )
