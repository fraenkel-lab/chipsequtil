#!/usr/bin/env python

import sys, re, os
from datetime import datetime
from optparse import OptionParser
from collections import defaultdict as dd
#from progressbar import ProgressBar
from csv import reader, writer

from chipsequtil import get_file_parts
from chipsequtil.util import MultiLineHelpFormatter as MF
from reStUtil import ReStDocument, ReStSimpleTable

usage = "%prog [options] <filename> [<filename>...]"
description="""\
Outputs various stats about the GERALD formatted file(s) input. If multiple
files are provided statistics are aggregated according to the specified output
format.  Output formats available via --format=X :

  # *python* - print an eval()'able python dictionary w/ counts
  # *rst* - print statistics in a reStructured text table (default)
  # *tab* - print statistics in a tab delimited form w/ header names

Except for *python* format, each input file has its own output line.  *python*
summarizes all alignments.
"""

parser = OptionParser(usage=usage,description=description,formatter=MF())
parser.add_option('--output',dest='output',default=None,help='write output to file [default: stdout]')
parser.add_option('--format',dest='format',type='choice',choices=['python','rst','tab'],default='rst',help='format to print out stats [default: %default]')

def log(st) :
    print datetime.now().isoformat()+' - '+st

re_digits_nondigits = re.compile(r'\d+|\D+')
def format_with_commas(value,format='%s'):
    parts = re_digits_nondigits.findall(format % (value,))
    for i in xrange(len(parts)):
        s = parts[i]
        if s.isdigit():
            parts[i] = _commafy(s)
            break
    return ''.join(parts)

def _commafy(s):

    r = []
    for i, c in enumerate(reversed(s)):
        if i and (not (i % 3)):
            r.insert(0, ',')
        r.insert(0, c)
    return ''.join(r)

if __name__ == '__main__' :

    opts,args = parser.parse_args(sys.argv[1:])

    gerald_fns = args

    all_stats = dd(int)
    stat_dicts = {}
    stats_fields = ["sample",
                    "total alignments",
                    "% align unique",
                    "# reads aligned unique",
                    "% align repeat",
                    "# reads align repeat",
                    "% align none",
                    "# reads align none"
                   ]


    data_rows = []
    for gerald_fn in gerald_fns :
        stats = stat_dicts[gerald_fn] = dd(int)

        fnpath,fn,fnbase,fnext = get_file_parts(gerald_fn)
        gerald_lines = reader(open(gerald_fn),delimiter='\t')
        for row in gerald_lines :
            m = re.match('^(\d+):(\d+):(\d+)$',row[10])
            if m is not None :
                stats['multiread'] += 1
                all_stats['multiread'] += 1
            else :
                stats[row[10]] += 1
                all_stats[row[10]] += 1

        tot_reads = sum(stats.values())/1.-stats.get('QC',0)
        unique_reads = sum([v for k,v in stats.items() if k.startswith('chr')])
        repeat_reads = stats.get('multiread',0)
        nomap_reads = stats.get('NM',0)
        data_row = [fn,format_with_commas(int(tot_reads)),
                    '%.1f'%(unique_reads/tot_reads*100),format_with_commas(unique_reads),
                    '%.1f'%(repeat_reads/tot_reads*100),format_with_commas(repeat_reads),
                    '%.1f'%(nomap_reads/tot_reads*100),format_with_commas(nomap_reads)]

        data_rows.append(data_row)

    out_f = open(opts.output,'w') if opts.output is not None else sys.stdout

    if opts.format == 'python' :
        out_f.write(dict(all_stats))
    elif opts.format == 'rst' :
        doc = ReStDocument(out_f)
        table = ReStSimpleTable(header=stats_fields,data=data_rows)
        doc.add(table)
        doc.write()
    elif opts.format == 'tab' :
        out_w = writer(out_f,delimiter='\t')
        out_w.writerow(stats_fields)
        out_w.writerows(data_rows)

    if opts.output is not None : out_f.close()
