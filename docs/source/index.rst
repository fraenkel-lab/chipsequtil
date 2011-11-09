.. ChIPSeqUtil documentation master file, created by
   sphinx-quickstart on Mon Oct 31 13:12:52 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ChIPSeqUtil's documentation!
=======================================

ChIPSeqUtil is a python module and accompanying set of scripts used in the
analysis of ChIPSeq short read data.  It is designed as a 'push-button' solution
that is easy for non-linux-experts to use but is flexible and extensible enough
to accomodate special cases when they inevitably arise. The default pipeline
performs the following analysis steps:

1. runs a peak caller (MACS by default)
2. optionally creates and stages bigwig files for viewing on UCSC Genome Browser
3. filters peaks based on confidence criteria (e.g. p-value)
4. maps peaks to genes using UCSC knownGene annotations
5. performs hypothesis-based motif analysis using TRANSFAC motifs
6. builds a web page consolidating results

ChIPSeqUtil has the following dependencies:

  - MACS (or some other peaks caller)
  - TAMO
  - reStUtil
  - pypeline
  - bx python

.. note:: add links to these bullets

ChIPSeqUtil has only been tested on ubuntu-based linux distributions and no
certification is made for other OSes.  That being said, some/all of it may
still work.

Contents:

.. toctree::
   :maxdepth: 2

   quick_start
   script_reference
   module_reference

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

