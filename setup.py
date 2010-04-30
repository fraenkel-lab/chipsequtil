#!/usr/bin/env python
from distutils.core import setup

setup(name='chipseq',
		version='0.1',
		scripts= ['bin/chipseq_pipeline.py',
		          'bin/combine_gerald_stats.py',
		          'bin/gerald_stats.py',
		          'bin/gerald_to_bed.py',
		          'bin/split_file.py',
		          'bin/split_qsub.py',
		          'bin/wait_for_qsub.py']
		)
