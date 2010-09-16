#!/usr/bin/env python

import re
import sys

from csv import reader, writer
from collections import defaultdict as dd
from optparse import OptionParser

from chipsequtil.util import MultiLineHelpFormatter as MF

usage = '%prog [options] <mapped known genes file>'
description = """Filter columns and rows from *join_mapped_known_genes.py* output which was \
invoked with *--binary-plus* and *--field-types* flags.  Specify full column names for either \
binding or expression data with the *--bind-cols* and *--affy-cols* arguments, respectively. \
The special fieldname *MAPPED* from *join_mapped_known_genes.py* is used to determine whether \
a file contains a mapping for each gene.  To filter genes by their associated binding or \
expression data, specify *--bind-filter* or *--affy-filter* as follows:

  - *any* - report gene if at least one input file maps to the gene
  - *all* - report gene if every input file maps to the gene
  - *absent* - report gene if no input file maps to the gene
  - *none* - do not filter genes at all (default)

Results of binding and expression filters are 'and'ed together, e.g.:

--bind-filter=all --affy-filter=absent

returns only genes for which all binding files and none of the expression files map.
"""
epilog='Note: when specifying column names, be sure to escape characters like (,),&,*,etc... \
that shells interpret with a \\, e.g. --bind-cols=-10\\*log10\\(pvalue\\)'
parser = OptionParser(usage=usage,description=description,epilog=epilog, formatter=MF())
parser.add_option('--bind-cols',dest='bind_cols',default='',help='comma delimited list of binding data column names to include, [default: all]')
parser.add_option('--affy-cols',dest='affy_cols',default='',help='comma delimited list of expression data column names to include, [default: all]')
parser.add_option('--bind-filter',dest='bind_filt',type='choice',choices=['any','all','absent','none'],default='none',help='gene set to include based on binding data [default: %default]')
parser.add_option('--affy-filter',dest='affy_filt',type='choice',choices=['any','all','absent','none'],default='none',help='gene set to include based on expression data [default: %default]')
parser.add_option('--output',dest='output',default=None,help='write output to file')


def match_headers(patts,field) :
    for p in patts :
        if field.endswith(p) : return True
    return False

def filter_vector(type,vec) :
    if type == 'any' :
        return '1' in vec
    elif type == 'all' :
        return all([x=='1' for x in vec])
    elif type == 'absent' :
        return not ('1' in vec)
    else :
        return True

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 1 :
        parser.error('Exactly one mapped file must be provided')

    map_fn = args[0]

    map_reader = reader(open(map_fn),delimiter='\t')
    headers = map_reader.next()
    bind_headers = [i for i,x in enumerate(headers) if x.startswith('BIND:')]
    bind_map_headers = [i for i,x in enumerate(headers) if x.startswith('BIND:') and x.endswith('MAPPED')]
    affy_headers = [i for i,x in enumerate(headers) if x.startswith('AFFY:')]
    affy_map_headers = [i for i,x in enumerate(headers) if x.startswith('AFFY:') and x.endswith('MAPPED')]

    if len(bind_headers) == 0 and len(affy_headers) == 0 :
        parser.error('No BIND: or AFFY: columns were found in the mapping, was *join_mapped_known_genes.py* run with the *--field-types* option?')

    # figure out which columns user wants
    header_indices = [0,1] # always output knowngene and symbol

    bind_header_patts = opts.bind_cols.split(',')
    header_indices += [i for i in bind_headers if match_headers(bind_header_patts,headers[i])]

    affy_header_patts = opts.affy_cols.split(',')
    header_indices += [i for i in affy_headers if match_headers(affy_header_patts,headers[i])]

    out_f = open(opts.output,'w') if opts.output else sys.stdout
    map_writer = writer(out_f,delimiter='\t')

    map_writer.writerow([headers[i] for i in header_indices])
    for rec in map_reader :
        bind_vector = [rec[i] for i in bind_map_headers]
        bind_pass = filter_vector(opts.bind_filt,bind_vector)

        affy_vector = [rec[i] for i in affy_map_headers]
        affy_pass = filter_vector(opts.affy_filt,affy_vector)

        if bind_pass and affy_pass :
            map_writer.writerow([rec[i] for i in header_indices])

    if opts.output : out_f.close()
