#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from ConfigParser import ConfigParser

usage = '%prog [options] <org key>'
description='Tool for retrieving sets of organism-specific settings and paths. Allows output of settings in a variety of shell environment syntaxes.'

if __name__ == '__main__' :

  opts, args = parser.parse_args(sys.argv[1:])
