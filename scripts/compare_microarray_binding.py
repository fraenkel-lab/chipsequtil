#!/usr/bin/env python

import sys

from csv import reader, writer
from collections import defaultdict as dd
from optparse import OptionParser
from subprocess import Popen, PIPE

from chipsequtil import MACSOutput, BEDOutput, AffyBiocFile

usage = '%prog -m <mapped MACS peaks file>|-b <mapped BED peaks file>|-a <mapped microarray file> [-m <MACS peaks file> ...] [-b <mapped BED peaks file> ...] [-a <mapped microarray file> ...]'
description = """Join all files on the first column, concatenating records with \
matching entries onto one line per entry.  Understands MACS peaks data as mapped \
with *map_peaks_to_known_genes.py* utility microarray data as mapped by \
*probeset_to_known_genes.py* utility, passed to program using *-m* and *-a* options \
respectively. Output is a file where genes with binding data (MACS, BED files) have \
column with a 1, 0 otherwise, and genes with microarray expression values have logFC \
and adjusted p-value colums for each microarray file input. Internally, uses \
*join_mapped_known_genes.py* with --binary-plus option to perform mapping and parses \
output.  MACS fields are listed first, followed by BED fields, followed by microarray \
fields."""

epilog="Note: microarray files should have been created by bioconductor, and all files should have a row of fieldnames as the first line"
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('-a','--affy-file',dest='affy_file',action='append',default=[],help='add a mapped microarray file')
parser.add_option('-b','--bed-file',dest='bed_file',action='append',default=[],help='add a mapped BED formatted peaks (*.bed) file')
parser.add_option('-m','--macs-file',dest='macs_file',action='append',default=[],help='add a mapped default MACS formatted peaks (*.xls) file')
parser.add_option('--output',dest='output',default=None,help='file to output joined records to [default: stdout]')

if __name__ == '__main__' :

    opts,args = parser.parse_args(sys.argv[1:])

    if len(args) > 0 :
        parser.error('There were non-option command line arguments passed, all files should have a preceeding option indicating filetype')

    if len(opts.macs_file) == 0 and len(opts.affy_file) == 0 :
        parser.error('No files were passed in, aborting')

    # call join_mapped_known_genes.py
    fn_map = {}
    fn_map['macs'] = ' '.join(['-m %s'%fn for fn in opts.macs_file])
    fn_map['bed'] = ' '.join(['-b %s'%fn for fn in opts.bed_file])
    fn_map['array'] = ' '.join(['-a %s'%fn for fn in opts.affy_file])
    join_call = 'join_mapped_known_genes.py --binary-plus %(macs)s %(bed)s %(array)s'%fn_map
    p = Popen(join_call, shell=True, stdout=PIPE,stderr=PIPE)
    stdout, stderr = p.communicate()
    if len(stderr) != 0 :
        print stderr

    joined_output = stdout.split('\n')
    joined_output = joined_output[:-1] if joined_output[-1] == '' else joined_output

    # determine which fields will end up in the file
    header = joined_output[0].split('\t')

    # always want gene and symbol
    field_indices = [0,1]

    # macs and bed fields are named by filename
    for fn in opts.macs_file+opts.bed_file :
        field_indices.append(header.index(fn))

    # affy fields are index(fn)+5, index(fn)+8
    for fn in opts.affy_file :
        # just add all the microarray columns
        fn_header_indices = [i for i,x in enumerate(header) if x.find(fn) != -1]
        field_indices.extend(fn_header_indices)

        #field_indices.append(header.index(fn))
        #field_indices.append(header.index(fn)+5)
        #field_indices.append(header.index(fn)+8)

    out_f = open(opts.output,'w') if opts.output else sys.stdout
    for line in joined_output :
        line = line.split('\t')
        out_f.write('\t'.join([line[i] for i in field_indices])+'\n')

    if opts.output :
        out_f.close()
