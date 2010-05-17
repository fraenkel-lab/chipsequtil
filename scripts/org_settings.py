#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from ConfigParser import ConfigParser, NoSectionError
from pprint import pformat
from chipsequtil import get_org_settings, get_all_settings

usage = '%prog [options] [<org key> [<org setting>]]'
description='Tool for retrieving sets of organism-specific settings and paths. Original paths are set at install time, and can be overridden in the file ~/.org_settings.cfg. Allows (eventually, not done yet) output of settings in a variety of shell environment syntaxes.  When run with an argument, settings are written to stdout by default in python pprint.pprint() repr() format.  When run without an argument, returns a listing of all settings available.'
parser = OptionParser(usage=usage,description=description)

if __name__ == '__main__' :

  opts, args = parser.parse_args(sys.argv[1:])

  # output depends on number of arguments passed
  output = ''

  # return everything we know about
  if len(args) == 0 :
    settings = get_all_settings()
    output = pformat(settings)

  # return all records from the specific organism
  elif len(args) in (1,2) :

    # make sure our config files have the requested organism
    try :
      settings = get_org_settings(args[0])
    except NoSectionError :
      sys.stderr.write('No entry %s found, available:\n'%args[0]+pformat(get_all_settings().keys())+'\nExiting\n')
      sys.exit(1)

    # return the requested field from the specific organism
    if len(args) == 2 :

      # make sure the config file has the setting for this organism
      try :
        output = pformat(settings[args[1]])
      except KeyError :
        sys.stderr.write('Setting %s not found for %s, choices:\n'+pformat(settings.keys())+'\nExiting\n')
        sys.exit(2)
    else :
      output = pformat(settings)
  else :
    parser.error('Provide zero, one, or two argments, found %s'%args)

  # bon voyage
  sys.stdout.write(output+'\n')

