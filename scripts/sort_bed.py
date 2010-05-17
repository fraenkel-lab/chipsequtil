#!/usr/bin/env python
import sys, os
from optparse import OptionParser
from collections import defaultdict as dd
from csv import reader, writer


usage = "%prog [options] <BED file> [<BED file> <BED file>...]"
description = """\
Sort the BED formatted files first by chromosome (field 1) and then by start
coordinate (field 2).  Lines from all files submitted are concatenated and
sorted in the final output."""
parser = OptionParser(usage=usage,description=description)
parser.add_option('--output',dest='output',default=sys.stdout,help='filename to write the sorted BED lines [default: stdout]')

if __name__ == '__main__' :

	opts, args = parser.parse_args(sys.argv[1:])

	if len(args) == 0 :
		parser.error("Must provide at least one file")

	fns = args
	chromos = dd(list)

	# load each chromosome separately
	for fn in fns :
		bed_reader = reader(open(fn),delimiter='\t')
		for line in bed_reader :
			chromos[line[0]].append(line)

	# determine where we're writing to
	if opts.output != sys.stdout :
		f = open(opts.output,'w')
	else :
		f = opts.output

	# write the chromos in lexicographic sorted order
	bed_writer = writer(f,delimiter='\t')
	for k in sorted(chromos.keys()) :

		# sort each chromosome's BED lines by stat position
		chromos[k].sort(key=lambda x: int(line[1]))
		bed_writer.writerows(chromos[k])
