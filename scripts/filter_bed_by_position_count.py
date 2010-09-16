#!/usr/bin/env python

import sys

from csv import reader, writer
from optparse import OptionParser

usage = '%prog [options] <bed file>'
description = """Analyze BED file and filter out alignments above some threshold \
that align to a single genomic position."""
epilog="Note: only works if BED file is sorted!"
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('-n','--max-count',dest='max_count',default=5,type='int',help='max tag count at a given position, filter above [default: %default]')
parser.add_option('--output',dest='output',default=None,help='write output to file')

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 1 :
        parser.error('Exactly one sorted .bed file is required')

    bed_fn = args[0]

    bed_reader = reader(open(bed_fn),delimiter='\t')
    out_f = open(opts.output,'w') if opts.output else sys.stdout
    bed_writer = writer(out_f,delimiter='\t')

    curr_key, curr_key_count = None, 0
    for rec in bed_reader :
        key = rec[:3] # chromosome, start, end
        if key != curr_key :
            curr_key, curr_key_count = key, 0
        if curr_key_count < opts.max_count :
            bed_writer.writerow(rec)
            curr_key_count += 1
        else :
            continue
    if opts.output : out_f.close()
