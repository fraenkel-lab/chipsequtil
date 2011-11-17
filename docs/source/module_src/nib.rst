
.. module:: chipsequtil.nib

nibFrag API
===========

These functions and classes are a native implementation of Jim Kent's nibFrag
utility and file format.  The scripts and classes read *.nib* files and can
extract sequences from them as fast or faster than the standalone tools, and
also make sequence data accessible and efficient from within python scripts.
There is no provided utility to create *.nib* files, the original source scripts
must be used and are not provided in this distribution.  They might be found on
`Jim Kent's homepage <http://users.soe.ucsc.edu/~kent/>`_.


The NibDB Class
---------------

.. autoclass:: NibDB
    :members:

Functions
---------

Most of these functions should not be used directly, rather they are called
by the NibDB class and implement the gritty details of reading *.nib* files.
Use the NibDB class instead unless you know what you're doing.


.. autofunction:: get_nib
.. autofunction:: get_nib_batch
.. autofunction:: get_nib_seq
.. autofunction:: get_nib_header
.. autofunction:: get_nib_header_batch
.. autofunction:: validate_nib_file
.. autofunction:: get_nib_seq_batch
