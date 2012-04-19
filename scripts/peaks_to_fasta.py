#!/usr/bin/env python

import os
import sys
import textwrap
import warnings
from optparse import OptionParser

from chipsequtil import BEDFile, MACSFile, get_file_parts, get_org_settings
from chipsequtil.nib import NibDB
from chipsequtil.sampling import rejection_sample_bg
from chipsequtil.util import MultiLineHelpFormatter
from chipsequtil.seq import write_fasta_to_file


usage='%prog [options] <organism> <peak file> [<peak file> ...]'
description='''Extract sequences for peaks in provided peak file(s).  Can \
interpret MACS or BED output, determined automatically by .xls or .bed extensions \
respectively (force explicit format with --peak-format option).  Outputs fasta \
sequences for the peaks in all files extracted from the reference genome specified \
by the output of *org_settings.py <organism> genome_dir* to stdout by default.\
Chromosome names in peak files must match nib filenames without extension (e.g. \
peak line: chr1 0  100 searches *genome_dir*/chr1.nib).  Fasta records have the \
following format:

><chromosome>:<start>-<end>;fn=<name of file>:<line number>;db_fn=<db filename>;fmt=<format>;<source alignment info>
<sequence...>

<db filename> is the filename where the sequence was extracted, <format> is the \
format of the input file (MACS or BED), and <source alignment info> contains all \
the fields from the originating alignment according to the source format.'''
parser = OptionParser(usage=usage,description=description,formatter=MultiLineHelpFormatter())
parser.add_option('--min-header',dest='min_header',action='store_true',help='only store <chromosome>:<start>-<end> in header')
parser.add_option('--peak-format',dest='peak_format',type='choice',
                  choices=['auto','MACS','BED'],default='auto',
                  help='peak file format, \'auto\' determines format by extension, choices: MACS, BED, auto [default: %default]')
parser.add_option('--output',dest='output',default=None,help='filename to output fasta records to [default: stdout]')
parser.add_option('--fixed-peak-width',dest='fixed_peak_width',type='int',default=None,help='return a fixed number of bases flanking peak summit (*summit* field in MACS, (end-start)/2 in BED), ignoring start/stop coords [default: None]')
parser.add_option('--wrap-width',dest='wrap_width',type='int',default=70,help='wrap fasta sequences to specified width. -1 indicates no wrap [default: %default]')


def bed_to_fasta(fn,db,min_header=False) :
    #headers,seqs = db.get_fasta_from_bed(fn)
    fastas = []
    bed_recs = BEDFile(fn)
    for i,rec in enumerate(bed_recs) :

        if opts.fixed_peak_width :
            midpoint = (rec['chromEnd']+rec['chromStart'])/2
            start = max(0,midpoint-opts.fixed_peak_width/2)
            end = min(midpoint+opts.fixed_peak_width/2,db.db_info[rec['chrom']]['nbases'])
            coords = start, end
        else :
            coords = start,end = int(rec['chromStart']), int(rec['chromEnd'])

        seq = db.get_seq(rec['chrom'], start, end)
        seq_fn = db.db_info[rec['chrom']]['path']

        header = '%s:%s;'%(rec['chrom'],'%d-%d'%(start,end))
        if not min_header :
            header = header.strip()+'%s:%d;fmt=BED;'%(fn,i)+ \
                     ';'.join(['%s=%s'%(k,str(v)) for k,v in rec.items()])
        fastas.append((header,seq))

    return fastas


def macs_to_fasta(fn,db,min_header=False) :
    macs_recs = MACSFile(fn)
    fasta = []
    for i,rec in enumerate(macs_recs) :

        if opts.fixed_peak_width :
            # adjust start and end peak position based on summit, ensuring we don't step outside of the reference sequence bounds
            start = max(0, rec['start']+rec['summit']-opts.fixed_peak_width/2)
            end = min(rec['start']+rec['summit']+opts.fixed_peak_width/2, db.db_info[rec['chr']]['nbases'])
            coords = start, end
        else :
            start, end = coords = rec['start'], rec['end']

        seq = db.get_seq(rec['chr'],start,end)
        seq_fn = db.db_info[rec['chr']]['path']

        header = '%s:%s'%(rec['chr'],'%d-%d'%coords)
        if not min_header :
            header += ';%s:%d;db_fn=%s;fmt=MACS;'%(fn,i,seq_fn) + \
                     ';'.join(['%s=%s'%(k,str(v)) for k,v in rec.items()])
        fasta.append((header,seq))

    return fasta


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 2 :
        parser.error('Must provide at least two non-option arguments')

    # instantiate the NibDB from the provided directory
    organism = args[0]
    nib_dir = get_org_settings(organism)['genome_dir']
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
        if opts.wrap_width == -1 :
            opts.wrap_width = sys.maxint
        write_fasta_to_file(dict(fasta_recs),opts.output,linelen=opts.wrap_width)
    else :
        for header, seq in fasta_recs :
            if opts.wrap_width != -1 :
                seq = textwrap.fill(seq,opts.wrap_width)
            sys.stdout.write('>%s\n%s\n'%(header,seq))
