#!/usr/bin/env python

import sys, re, os
from optparse import OptionParser
from collections import defaultdict as dd

parser = OptionParser()

if __name__ == '__main__' :

	opts, args = parser.parse_args(sys.argv[1:])

	all_stats = dd(int)
	for fn in args :
		d = eval(open(fn).read())
		for k,v in d.items() :
			all_stats[k] += v
			all_stats['tot. aligns'] += v

	keys = all_stats.keys()
	keys.sort()
	keys.remove('tot. aligns')

	for k in keys :
		print k,':',all_stats[k],'(%.2f)'%(float(all_stats[k])/all_stats['tot. aligns'])

	print 'tot. aligns',':',all_stats['tot. aligns']
