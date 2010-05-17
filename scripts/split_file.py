#!/usr/bin/env python

from optparse import OptionParser
from datetime import datetime
from subprocess import Popen, PIPE
import itertools
import sys, os, getpass, re

usage = "[%prog] [options] filename"
description = """\
Split <filename> into a set of files with either a specific number of lines
(--split-type=lines, default) or into a specific number of files (--split-type=
count).  Files are created with .XXXX appended, indicating the number of file
split. Writes files to current working directory unless otherwise specified.
"""

parser = OptionParser(usage=usage,description=description)
parser.add_option('--type',dest='split_type',type='choice',choices=['lines','count'],default='lines',help='how to split the file (WARNING: count does not preserve the sequence of lines in the original file when splitting) [default: %default]')
#parser.add_option('--split-arg',dest='split_arg',default='1000',help='integer argument for split type (size specified as Xb, XK, XM, or XG, others are integers) [default: %default]')
parser.add_option('--arg',dest='split_arg',type='int',default=1000,help='integer argument for split type [default: %default]')
parser.add_option('--outdir',dest='outdir',default='.',help='directory to put the split files in [default: %default]')

def get_file_parts(fn) :
    fpath,fname = os.path.split(fn)
    fbase,fext = os.path.splitext(fname)
    return fpath,fname,fbase,fext

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 1 :
        parser.print_usage()
        sys.exit(1)

    filename = args[0]
    abs_filename = os.path.abspath(filename)

    # check to ensure filename exists
    if not os.path.exists(abs_filename) :
        sys.stderr.write('File %s does not exist, exiting\n'%abs_filename)
        parser.print_usage()
        sys.exit(2)

    # split the file
    split_size = opts.split_arg
    fpath,fname,fbase,fext = get_file_parts(abs_filename)
    if opts.split_type == 'lines' :
        curr_split = 0 # for first condition
        split_fd = None
        for i,l in enumerate(open(abs_filename)) :
            if i%split_size == 0 :
                if split_fd : split_fd.close() # close it if we aren't on the first split
                split_fd = open(os.path.join(opts.outdir,fname)+'.%04d'%curr_split,'w')
                curr_split += 1
            split_fd.write(l)
        nlines = i
    elif opts.split_type == 'count' :
        # create split_size split files by writing lines round robin
        split_fds = [open(os.path.join(opts.outdir,fname)+'.%04d'%x,'w') for x in range(split_size)]
        split_cycle = itertools.cycle(split_fds)
        for i,l in enumerate(open(abs_filename)) :
            split_cycle.next().write(l)
        nlines = i

        # close all the handles
        [fd.close() for fd in split_fds]

    elif opts.split_type == 'size' :
        # parse split_arg argument, into integer if split_type is 'size'
        if opts.split_type == 'size' :
            m = re.match('^(\d+)([bKMG])$',opts.split_arg)
            if m is None :
                sys.stderr.write("Incorrect --split-arg argument for --split-type=size, I understand only X[bKMG], exiting\n")
                parser.print_usage()
                sys.exit(3)
            else :
                size_d = {'b':1,'K':1024,'M':pow(1024,2),'G':pow(1024,3)}
                split_size = int(m.groups()[0])*size_d[m.groups()[1]]

        fd = open(abs_filename)
        curr_split_size = 0

