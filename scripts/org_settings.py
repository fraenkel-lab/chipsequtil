#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
from ConfigParser import ConfigParser, NoSectionError
from pprint import pformat

from chipsequtil import get_org_settings, get_global_settings, get_all_settings, get_local_settings, GLOBAL_SETTINGS_FN, LOCAL_SETTINGS_FN

usage = '%prog [options] [<org key> [<org setting>]]'
description='Tool for retrieving sets of organism-specific settings and paths. \
Original paths are set at install time, and can be overridden in the file ~/.org_\
settings.cfg. Allows (eventually, not done yet) output of settings in a variety \
of shell environment syntaxes.  When run with an argument, settings are written \
to stdout by default in python pprint.pprint() repr() format.  When run without \
an argument, returns a listing of all settings available.'
parser = OptionParser(usage=usage,description=description)
parser.add_option('-s','--syntax',dest='syntax',type='choice',\
                  choices=['python','bash'],default='python',help='syntax flavor \
                  of output to produce [default: %default]')
parser.add_option('-l','--list',dest='list_sets',action='store_true',help='print \
                  all available settings for human consumption')


def obj_to_format(obj,format='python') :
    '''Convert *obj* into a string that can be evaluated in the environment \
    indicated in *format*.

    obj -- a string, a dict of values, or a dict of dicts of values
    format -- python (default), or bash
    '''

    r = ''
    if format == 'python' :
        r = pformat(obj)
    elif format == 'bash' :
        statements = []
        export_tmpl = 'export %s=%s'

        # should only get a string, a dict
        if isinstance(obj,str) :
            return obj
        # dict
        elif isinstance(obj,dict) :
            for k1, v1 in obj.items() :
                # dict of dicts
                if isinstance(v1,dict) :
                    # these should be literal values
                    for k2, v2 in v1.items() :
                        statements.append(export_tmpl%('_'.join([k1,k2]).upper(),\
                                          str(v2)))
                elif isinstance(v1,str) :
                    s = '"'+v1+'"' if v1.count(' ') != 0 else str(v1)
                    statements.append(export_tmpl%(k1.upper(),str(s)))

        r = '\n'.join(statements)

    return r


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    # output depends on number of arguments passed
    output = ''

    # return everything we know about
    if len(args) == 0 :

        if opts.list_sets :
            # global settings
            settings = get_global_settings()
            output = 'Global settings: (%s)\n'%GLOBAL_SETTINGS_FN
            output += obj_to_format(settings,opts.syntax) + '\n'

            # local settings
            settings = get_local_settings()
            output += 'Local settings: (%s)\n'%LOCAL_SETTINGS_FN
            output += obj_to_format(settings,opts.syntax)
        else :
            settings = get_all_settings()
            output += obj_to_format(settings,opts.syntax)


    # return all records from the specific organism
    elif len(args) in (1,2) :

        # make sure our config files have the requested organism
        try :
            settings = get_org_settings(args[0])
        except NoSectionError :
            sys.stderr.write('No entry %s found, available:\n'%args[0]+\
                             pformat(get_all_settings().keys())+'\nExiting\n')
            sys.exit(1)

        # return the requested field from the specific organism
        if len(args) == 2 :

            # make sure the config file has the setting for this organism
            try :
                output = obj_to_format(settings[args[1]],opts.syntax)
            except KeyError :
                sys.stderr.write('Setting %s not found for %s, choices:\n'+
                                 pformat(settings.keys())+'\nExiting\n')
                sys.exit(2)
        else :
            output = obj_to_format(settings,opts.syntax)
    else :
        parser.error('Provide zero, one, or two argments, found %s'%args)

    # bon voyage
    sys.stdout.write(output+'\n')

