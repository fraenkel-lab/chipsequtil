#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
from ConfigParser import ConfigParser, NoSectionError
from pprint import pformat

from chipsequtil import get_org_settings, get_global_settings, get_all_settings, get_local_settings, GLOBAL_SETTINGS_FN, LOCAL_SETTINGS_FN

usage = '%prog [options] [<org key> [<org setting>]]'
description='''Tool for retrieving sets of organism-specific settings and paths.
Original paths are set at install time, and can be overridden in the file ~/.org
settings.cfg. Allows output of settings in a variety of shell environment
syntaxes.  The tool attempts to guess which shell environment is being used by
examining the SHELL environment variable unless explicitly set.  When run without
an argument, returns a listing of all settings available.
'''
parser = OptionParser(usage=usage,description=description)
parser.add_option('-s','--syntax',dest='syntax',type='choice',\
                  choices=['auto','python','bash','tcsh'],default='auto',help='syntax flavor \
                  of output to produce [default: %auto]')
parser.add_option('-l','--list',dest='list_sets',action='store_true',help='print \
                  all available settings for human consumption')


def obj_to_format(obj,format='python') :
    '''Convert *obj* into a string that can be evaluated in the environment \
    indicated in *format*.

    obj -- a string, a dict of values, or a dict of dicts of values
    format -- python (default), or bash
    '''

    if format == 'auto' :
        format = os.environ.get('SHELL','python').split('/')[-1]

    r = ''
    if format == 'python' :
        r = pformat(obj)
    elif format in ['sh','bash','zsh','csh','tcsh'] :
        statements = []
        if format in ['sh','bash','zsh'] :
            export_tmpl = 'export %s=%s'
        elif format in ['csh','tcsh'] :
            export_tmpl = 'setenv %s %s'

        # dict
        if isinstance(obj,dict) :
            for k1, v1 in obj.items() :
                # dict of dicts
                if isinstance(v1,dict) :
                    # these should be literal values
                    for k2, v2 in v1.items() :
                        statements.append(export_tmpl%('_'.join([k1,k2]).upper(),\
                                          str(v2)))
                else :
                    v1 = str(v1)
                    s = '\''+v1+'\'' if v1.count(' ') != 0 else str(v1)
                    statements.append(export_tmpl%(k1.upper(),str(s)))
        else :
            return str(obj)

        r = '\n'.join(statements)

    return r


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    # output depends on number of arguments passed
    output = ''

    # return everything we know about
    if len(args) == 0 :

        if opts.list_sets :

            # always use python formatting when listing
            opts.syntax = 'python'

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
                sys.stderr.write('Setting %s not found for %s, choices:\n'%(args[1],args[0])+
                                 pformat(settings.keys())+'\nExiting\n')
                sys.exit(2)
        else :
            output = obj_to_format(settings,opts.syntax)
    else :
        parser.error('Provide zero, one, or two argments, found %s'%args)

    # bon voyage
    sys.stdout.write(output+'\n')

