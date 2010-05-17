#!/usr/bin/env python
import sys, os
from distutils.core import setup, run_setup
from distutils.command.build_py import build_py
from distutils.file_util import copy_file
from distutils.sysconfig import EXEC_PREFIX, PREFIX

# convenience is king
opj = os.path.join

# copy org_settings.cfg in current directory to source tree
org_settings_fn = 'org_settings.cfg'
org_settings_path = os.path.abspath(org_settings_fn)
dist_settings_path = opj(os.getcwd(),'src','chipsequtil',org_settings_fn)
if os.path.exists(os.path.abspath(org_settings_fn)) :
  copy_file(org_settings_path,dist_settings_path)
else :
  sys.stderr.write('WARNING: org_settings.cfg could not be found in distribution root directory.  org_settings.py script may not work properly.\n')
  open(dist_settings_path,'w').close() # touch the file so we don't get errors during install, user has been warned

scripts = ['scripts/chipseq_pipeline.py',
           'scripts/combine_gerald_stats.py',
           'scripts/gerald_stats.py',
           'scripts/gerald_to_bed.py',
           'scripts/sort_bed.py',
           'scripts/split_file.py',
           'scripts/split_qsub.py',
           'scripts/wait_for_qsub.py',
           'scripts/map_peaks_to_genes.py',
           'scripts/org_settings.py',
           'scripts/nibFrag.py']

# distutils doesn't handle uninstalling things, this class deletes all the files this package installs
# if it has appropriate permissions to do it, otherwise print out the files that must be deleted to uninstall
class uninstall(build_py) :
  def run(self) :

    #TODO consider doing this later, use the install -f|--force option for now

    # delete modules
    print self.distribution.py_modules

    # delete extensions
    print self.distribution.ext_modules

    # delete packages
    print self.distribution.packages

    # delete package data
    print self.distribution.package_data

    # delete scripts
    print self.distribution.scripts

    print self.distribution.get_command_obj('install').get_outputs()

  def remove_path(self,path) :
    '''Attempt to remove the specified path, returning non-zero status code on error'''

# setup and install
setup(name='chipsequtil',
      version='0.1',
      package_dir={'':'src'},
      py_modules=['chipsequtil.nib'],
      packages=['chipsequtil'],
      #package_data={'': ['org_settings.cfg']},
      scripts=scripts,
      cmdclass={'uninstall': uninstall},
     )
