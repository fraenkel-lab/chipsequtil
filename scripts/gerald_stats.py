#!/usr/bin/env python

import sys, re, os
from datetime import datetime
from optparse import OptionParser
from collections import defaultdict as dd
#from progressbar import ProgressBar
from csv import reader

usage = "%prog [options] <filename> [<filename>...]"
description="""\
Outputs various stats about the GERALD formatted file(s) input. If multiple
files are provided all statistics are summarized unless the --separate option
is supplied.
"""

parser = OptionParser(usage=usage,description=description)

def log(st) :
    print datetime.now().isoformat()+' - '+st

def get_file_parts(path) :
    path,fn = os.path.split(path)
    basename,ext = os.path.splitext(fn)
    return path,fn,basename,ext

if __name__ == '__main__' :

    opts,args = parser.parse_args(sys.argv[1:])

    gerald_fns = args
    #gerald_fns = ['/home/labadorf/sge/stat_test/xa']

    stats = dd(int)
    #log('Starting')
    for gerald_fn in gerald_fns :

        # hack to get the number of lines in the file for progressbar purposes...

        #log('Processing: '+gerald_fn)
        #nlines = max((i for i in enumerate(open(gerald_fn))))[0]+1
        #pbar = ProgressBar(maxval=nlines).start()
        gerald_lines = reader(open(gerald_fn),delimiter='\t')
        for row in gerald_lines :
            m = re.match('^(\d+):(\d+):(\d+)$',row[10])
            if m is not None :
                stats['multiread'] += 1
            else :
                stats[row[10]] += 1
            #pbar.update()
        #pbar.finish()

    print dict(stats)
