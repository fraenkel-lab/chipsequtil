
.. module:: chipsequtil.seq

Sequence data functions and classes
===================================

This module has simple methods for reading in FASTA and FASTQ formatted files.
*fasta_itr* and *fastq_itr* should be used when it is unnecessary or undesired
to have all sequences loaded into memory.  *FASTAFile* and *FASTQFile* classes
store all sequence information in memory, but allow efficient dictionary-style 
random access to sequences and quality scores as well as repeated whole-file
iteration.

Functions
---------

.. autofunction:: fasta_itr
.. autofunction:: fasta_to_dict
.. autofunction:: write_fasta_to_file

.. autofunction:: fastq_itr
.. autofunction:: fastq_to_dict
.. autofunction:: write_fastq_to_file

Classes
-------

.. autoclass:: FASTAFile
    :members:

.. autoclass:: FASTQFile
    :members:


