#!/usr/bin/env python
import sys, os
from optparse import OptionParser
from collections import defaultdict as dd


usage = "%prog [options] <BED file> [<BED file> <BED file>...]"
description = "\
Sort the BED formatted files first by chromosome (field 1) and then by start
coordinate (field 2).  Lines from all files submitted are concatenated and
sorted in the final output."
parser = OptionParser(usage=usage,description=description)

if __name__ == '__main__' :

	opts, args = parser.parse_args(sys.argv[1:])


