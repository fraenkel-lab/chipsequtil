#!/usr/bin/env python

import os
import re
import sys

from optparse import OptionParser
from csv import DictReader, DictWriter
from chipsequtil import get_file_parts, GERALDOutput

usage = "%prog [options] <GERALD file> [<GERALD file>...]"

description = """\
Convert the GERALD alignment formatted files into BED format.  Input file named
<path>/<filename>.<ext> is translated into <path>/<filename>.bed unless --output
or --stdout is specified, in which case formatted lines are written to file or
standard output, respectively.  If multiple input files are supplied with the
--output or --stdout option all formatted lines are concatenated together.
Formatting only occurs for GERALD input lines that have a valid Match Position
field (i.e. successfully aligned somewhere)."""

parser = OptionParser(usage=usage, description=description)
parser.add_option('--output',dest='output',default=None,help='write all records to file')
parser.add_option('--stdout',dest='stdout',action='store_true',help='write out all formatted lines to stdout')
parser.add_option('--min-fields',dest='min_fields',action='store_true',help='only format the first three fields')
parser.add_option('--pass-only',dest='pass_only',action='store_true',help='only format lines with Y in the Pass Filtering field')
parser.add_option('--chromo-strip',dest='chromo_strip',default='.fa',help='pattern to remove from chromo field in BED output (e.g. --chromo-strip=.fa to remve .fa from chrX.fa) [default: %default]')



if __name__ == '__main__' :

    opts,args = parser.parse_args(sys.argv[1:])

    if len(args) == 0 :
        parser.print_usage()
        sys.exit(1)

    gerald_fns = args

    # step through the files
    for gerald_fn in gerald_fns :
        path,fn,fnbase,fnext = get_file_parts(gerald_fn)
        bed_lines = []


        # where to write output to
        if opts.stdout :
            f_out = sys.stdout
        else :
            f_out = open(os.path.join(path,fnbase+'.bed'),'w')

        # process input
        gerald_d = DictReader(open(gerald_fn),fieldnames=GERALDOutput.FIELD_NAMES,delimiter='\t')
        for line_d in gerald_d :
            if (opts.pass_only and line_d['filtering'] == 'Y' and line_d['match_pos'] != '') or (not opts.pass_only and line_d['match_pos'] != '') :

                if opts.chromo_strip is not None :
                    line_d['match_chromo'] = line_d['match_chromo'].replace(opts.chromo_strip,'')

                outline = [line_d['match_chromo'], # chromosome
                           line_d['match_pos'], # start
                           str(int(line_d['match_pos'])+len(line_d['read'])), # end
                           line_d['read'], # read
                           '0', # score
                           '+' if line_d['match_strand'] == 'F' else '-', # strand
                           '-', # thickStart
                           '-', # thickEnd
                           '0,0,255' if line_d['match_strand'] == 'F' else '255,0,0', # itemRgb 
                          ]
                outline = '\t'.join(outline)
                f_out.write(outline+'\n')
                #bed_lines.append(bed)

        # this is the slow way
        #for line in open(gerld_fn) :
        #    grld = GERALDOutput(line)
        #    if (opts.pass_only and grld.filtering == 'Y' and grld.match_pos != '') or (not opts.pass_only and grld.match_pos != '') :
        #        bed = gerald_to_bed(grld,opts.min_fields)
        #        f_out.write(bed.output_format())
        #        #bed_lines.append(bed)

