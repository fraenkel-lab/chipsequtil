#!/usr/bin/env python

import sys

from collections import defaultdict
from csv import reader
from optparse import OptionParser

from bx.intervals.intersection import IntervalTree, Interval

usage = '%prog [options] <from> <to>'
description = """Find records in <to> interval file that map to records in
<from> interval file.  Files should be tab delimited and are expected to have
a chromosome column, a start column, and an end column.  The indices of these
columns can be specified on the command line but by default are the first
three columns, respectively.  Prints out to stdout by default one new line
separated row per row in <from> with a line from <to> where there is a mapping.
If no mapping is found (e.g. when specifying a maximum margin to search within)
the word None is printed.  By default only prints nearest record, with ties
settled by smallest line number in <to>."""
parser = OptionParser(usage=usage,description=description)
parser.add_option('-w','--window',dest='window',type="float",nargs=2,
                  default=(1e9,1e9),
                  help="window as <int upstream> <int downstream> to search for intervals [default: %default]")
parser.add_option('-f','--from',dest='from_ind',type="int",nargs=3,
                  default=(0,1,2),
                  help="coordinates of chromosome, start, stop in <from> file")
parser.add_option('-i','--skip-from-header',dest='skip_fh',action='store_true',
                  help="<from> has a header that should be skipped")
parser.add_option('-t','--to',dest='to_ind',type="int",nargs=3,
                  default=(0,1,2),
                  help="coordinates of chromosome, start, stop in <to> file")
parser.add_option('-j','--skip-to-header',dest='skip_th',action='store_true',
                  help="<to> has a header that should be skipped")

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 2 :
        parser.error('Exactly 2 non-option arguments are required')

    from_fn, to_fn = args

    chr_trees = defaultdict(IntervalTree)
    chr_sizes = defaultdict(lambda : dict(minstart=sys.maxint,maxend=0))

    if any([x > 1e9 for x in opts.window]) :
        parser.error('Window maximum is +/- 1e9')

    to_reader = reader(open(to_fn),delimiter='\t')
    if opts.skip_th :
        to_header = to_reader.next()

    to_chr, to_st, to_en = opts.to_ind
    for r in to_reader :
        i = Interval(int(r[to_st]),
                     int(r[to_en]),
                     value=r,
                     chrom=r[to_chr]
                     )
        chr_trees[r[to_chr]].insert_interval(i)
        chr_sizes[r[to_chr]]['minstart'] = min(int(r[to_st]),chr_sizes[r[to_chr]]['minstart'])
        chr_sizes[r[to_chr]]['maxend'] = max(int(r[to_st]),chr_sizes[r[to_chr]]['maxend'])

    # window default is 1e9 because no chromosome is more than
    # ten billion base pairs, right?!
    def find_nearest(t,s,e,window=(1e9,1e9)) :

        # look for record within intervals
        inside = t.find(s,e)
        
        if len(inside) >= 1 : # pick the first one, list returned is sorted
            return inside[0]

        i = Interval(s,e)
        before = t.upstream_of_interval(i,max_dist=window[0])
        after = t.downstream_of_interval(i,max_dist=window[1])

        before = before[0] if len(before) != 0 else None
        after = after[0] if len(after) != 0 else None

        if before and after :
            b_dist = min(abs(before.end-s),abs(e-before.start))
            a_dist = min(abs(after.end-s),abs(e-after.start))
            nearest = before if b_dist < a_dist else after
        elif before :
            nearest = before
        elif after :
            nearest = after
        else :
            nearest = None
        return nearest

    # now go through the from file
    from_reader = reader(open(from_fn),delimiter='\t')
    if opts.skip_fh : from_reader.next()

    from_chr, from_st, from_en = opts.from_ind
    if opts.skip_th :
        print '\t'.join(to_header)
    for r in from_reader :
        t = find_nearest(chr_trees[r[from_chr]],int(r[from_st]),int(r[from_en]),
                         window=opts.window)
        if t :
            print '\t'.join(t.value)
        else :
            print t
    """
    # tests
    print 'interval is before any other interval in tree'
    t = find_nearest(chr_trees['chr2'],10388500,10388510)
    print '\tCorrect answer: %s, Returned answer: %s'%('mmu-mir-466f-1',t.value),t
    print 'interval is after any other interval in tree'
    t = find_nearest(chr_trees['chr1'],200000000,200000010)
    print '\tCorrect answer: %s, Returned answer: %s'%('mmu-mir-29c',t.value),t
    print 'interval is between intervals'
    t = find_nearest(chr_trees['chr3'],89773941,89774021)
    print '\tCorrect answer: %s, Returned answer: %s'%('mmu-mir-190b',t.value),t
    print 'interval is inside another interval'
    t = find_nearest(chr_trees['chr3'],89873999,89874001)
    print '\tCorrect answer: %s, Returned answer: %s'%('mmu-mir-190b',t.value), t
    print 'interval is too far from anything to return anything'
    t = find_nearest(chr_trees['chr3'],89773941,89774021,window=10)
    print '\tCorrect answer: None, Returned answer: %s'%t
    """
