#!/usr/bin/env python

import os
import sys
import textwrap
import warnings
from optparse import OptionParser, IndentedHelpFormatter

from chipsequtil import BEDFile, MACSFile
from chipsequtil.nib import NibDB
from chipsequtil.sampling import rejection_sample_bg_chris_reimpl
from chipsequtil.util import MultiLineHelpFormatter
from TAMO.seq import Fasta


usage='%prog [options] <nib dir> <peak file> [<peak file> ...]'
description='''Extract sequences for peaks in provided peak file(s).  Can \
interpret MACS or BED output, determined automatically by .xls or .bed extensions \
respectively (force explicit format with --peak-format option).  Outputs fasta \
sequences for the peaks in all files, extracted from .nib files found in provided \
nib directory to stdout by default.  Chromosome names in peak files must match nib \
filenames without extension (e.g. peak line: chr1 0  100 searches \
<nib dir>/chr1.nib).  Fasta records have the following format:

><chromosome>:<start>-<end>;fn=<name of file>:<line number>;db_fn=<db filename>;fmt=<format>;<source alignment info>
<sequence...>

<db filename> is the filename where the sequence was extracted, <format> is the \
format of the input file (MACS or BED), and <source alignment info> contains all \
the fields from the originating alignment according to the source format.'''
parser = OptionParser(usage=usage,description=description,formatter=MultiLineHelpFormatter())
parser.add_option('--min-header',dest='min_header',action='store_true',help='only store <chromosome>:<start>-<end> in header')
parser.add_option('--peak-format',dest='peak_format',type='choice',
                  choices=['auto','MACS','BED'],default='auto',
                  help='peak file format, \'auto\' determines format by extension,\
                  choices: MACS, BED, auto [default: %default]')
parser.add_option('--output',dest='output',default=None,help='filename to output \
                  fasta records to [default: stdout]')
parser.add_option('--bg',dest='bg',default='rej_samp',type='choice',
                  choices=['rej_samp'],help='generate background sequences that match \
                  the peak sequence distribution')



def bed_to_fasta(fn,db,min_header=False) :
    bed_recs = BEDFile(fn)
    fasta = []
    for i,rec in enumerate(bed_recs) :
        seq = db.get_seq(rec['chrom'],rec['chromStart'],rec['chromEnd'])
        seq_fn = db.db_info[rec['chrom']]['path']
        header = '%s:%s-%s'%(rec['chrom'],rec['chromStart'],rec['chromEnd'])
        if not min_header :
            header += ';%s:%d;db_fn=%s;fmt=BED;'%(fn,i,seq_fn) +\
                      ';'.join(['%s=%s'%(k,str(v)) for k,v in rec.items()])
        fasta.append((header,seq))
    return fasta

def macs_to_fasta(fn,db,min_header=False) :
    macs_recs = MACSFile(fn)
    fasta = []
    for i,rec in enumerate(macs_recs) :
        seq = db.get_seq(rec['chr'],rec['start'],rec['end'])
        seq_fn = db.db_info[rec['chr']]['path']
        header = '%s:%s-%s'%(rec['chr'],rec['start'],rec['end'])
        if not min_header :
            header += ';%s:%d;db_fn=%s;fmt=MACS;'%(fn,i,seq_fn) +\
                     ';'.join(['%s=%s'%(k,str(v)) for k,v in rec.items()])
        fasta.append((header,seq))
    return fasta


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 2 :
        parser.error('Must provide at least two non-option arguments')

    # instantiate the NibDB from the provided directory
    nib_dir = args[0]
    nib_db = NibDB(nib_dirs=[nib_dir])

    # determine specified format
    peak_fmt = opts.peak_format

    peak_fns = args[1:]

    # determine if there is an output file
    if opts.output :
        out_f = open(opts.output,'w')
    else :
        out_f = sys.stdout

    fasta_recs = []
    for peak_fn in peak_fns :
        # if --peak-format is auto, figure format out from extension
        if opts.peak_format == 'auto' :
            fnbase, fnext = os.path.splitext(peak_fn)
            if fnext.lower() == '.bed' : # BED file
                peak_fmt = 'BED'
            elif fnext.lower() == '.xls' : # MACS file
                peak_fmt = 'MACS'
            else  :
                warnings.warn('Peak format specified as auto but file extension \
                               not recognized in file %s, skipping'%peak_fn)
                continue

        if peak_fmt == 'BED' :
            fasta_recs.extend(bed_to_fasta(peak_fn,nib_db,min_header=opts.min_header))
        elif peak_fmt == 'MACS' :
            fasta_recs.extend(macs_to_fasta(peak_fn,nib_db,min_header=opts.min_header))

    # write out foreground to file
    if opts.output :
        Fasta.write(dict(fasta_recs),opts.output)
    else :
        for header, seq in fasta_recs :
            sys.stdout.write('>%s\n%s\n'%(header,seq))

    # create bg sequences if requested
    if opts.bg :
        sys.stderr.write('running rejection sampling\n')
        bg_dict = rejection_sample_bg_chris_reimpl(dict(fasta_recs),'mouse')
        Fasta.write(bg_dict,'rej_samp_bg.fa')
