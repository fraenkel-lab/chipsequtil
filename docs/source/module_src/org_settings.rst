
The `org_settings` System
=========================

Many scripts in this package require a number of different source files that all
correspond to a single reference genome (*e.g.* mm9).  The `org_settings` set of
functions and *org_settings.py* script consolidates sets of paths/variables that
correspond to different references to be bundled together in a customizable,
accessible way.  The bundles are configured as a package-wide settings on install
and alternatively by a user-specific configuration file.  The format of the file
follows the conventions in `configparser`_.

.. _configparser: http://docs.python.org/library/configparser.html

Reference genomes are specified in a configuration file as follows::

    [mm9]
    description=UCSC mm9 (July '07 build) with full TRANSFAC hypothesis set
    genome=mm9
    genome_dir=/nfs/genomes/mouse_gp_jul_07
    genome_size=2107000000
    ucsc_chrom_sizes=%(genome_dir)s/%(genome)s.chrom.sizes
    annotation_path=%(genome_dir)s/anno/refFlat-%(genome)s.txt
    refgene_anno_path=%(genome_dir)s/anno/refFlat-%(genome)s.txt
    known_gene_anno_path=%(genome_dir)s/anno/knownGene-%(genome)s.txt
    known_gene_xref_path=%(genome_dir)s/anno/kgXref-%(genome)s.txt
    affy_to_known_path=%(genome_dir)s/anno/knownToMOE43-%(genome)s.txt
    theme_hypotheses=/nfs/vendata/cwng/TRANSFAC/2010_transfac_vert_all_filtic9.tamo
    theme_markov=/nfs/data/cwng/chipseq/hypotheses/Mouse.markov

This will make **mm9** available as an organism reference to the `org_settings`
functions. The *ucsc_chrom_sizes*, *annotation_path*, *refgene_anno_path*,
*known_gene_anno_path*, *known_gene_xref_path*, and *affy_to_known_path* are
files downloaded from http://hgdownload.cse.ucsc.edu/downloads.html organims
annotation databases.  The fields in the above example are all required for the
package to work properly - however, any additional variables may be added as
desired.

API Functions
-------------

.. module:: chipsequtil

.. autofunction:: get_org_settings
.. autofunction:: get_all_settings
.. autofunction:: get_global_settings
.. autofunction:: get_local_settings
.. autofunction:: check_org_settings

The *org_settings.py* script
----------------------------

The script *org_settings.py* is a command line interface into the `org_settings`
system.  It has the following usage::

  $> org_settings.py -h
  Usage: org_settings.py [options] [<org key> [<org setting>]]

  Tool for retrieving sets of organism-specific settings and paths. Original
  paths are set at install time, and can be overridden in the file ~/.org
  settings.cfg. Allows output of settings in a variety of shell environment
  syntaxes.  The tool attempts to guess which shell environment is being used by
  examining the SHELL environment variable unless explicitly set.  When run
  without an argument, returns a listing of all settings available.

  Options:
    -h, --help            show this help message and exit
    -s SYNTAX, --syntax=SYNTAX
                          syntax flavor                   of output to produce
                          [default: %auto]
    -l, --list            print                   all available settings for
                          human consumption
  $> org_settings.py -s bash mm9 genome_dir
  /nfs/genomes/mouse_gp_jul_07
  $>

If you use bash as your shell, you can use shell expansion to conveniently build
commands such as the following::

  $> map_peaks_to_known_genes.py $(org_settings.py mm9 known_gene_anno_path) \
     $(org_settings.py mm9 known_gene_xref_path) macs_peaks.xls

Installing
----------

The file *org_settings.cfg* exists in the root directory of the source distribution.
This file should be modified and then copied into the *src/chipsequtil/* directory
before installation for org settings that should be available on the system as a
whole.  Alternatively, users may create the file *.org_settings.cfg* in their home
directories and add sections like the one above so they may customize their own
sets of variables.
