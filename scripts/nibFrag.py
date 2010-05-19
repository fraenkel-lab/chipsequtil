#!/usr/bin/env python
# nibFrag.py - a python implementation of Jim Kent's nibFrag command line utility

import sys
import warnings
from optparse import OptionParser, OptionGroup

from chipsequtil import get_file_parts, BEDFile
from chipsequtil.nib import get_nib_seq_batch, validate_nib_file, NibException, NOMASK, MASK, HARDMASK

usage = '%prog [options] file.nib start end strand [outfile]\n  -- or --\n%prog [options] --batch file.nib batchfile [batchfile ...]'
description = """A python implementation of Jim Kent's nibFrag utility that allows outputting to \
stdout.  Otherwise the functionality is identical for the non-batch usage.  Batch mode accepts \
one or more files containing sets of coordinates to extract from the nib file.  Only BED formatting \
is accepted at the moment. All sequences are concatenated together in FASTA format.  To retrieve the \
entire sequence, use END as the end argument."""
epilog="Note: When specifying --name optionin batch mode, also specify --dbHeader to ensure unique FASTA headers."
parser = OptionParser(usage=usage,description=description,epilog=epilog)
#parser.add_option('--output',dest='output',default=sys.stdout,help='filename to write output to [default: stdout]')
parser.add_option('--no-header',dest='no_header',action='store_true',help='only output sequence (no fasta header)')
parser.add_option('--wrap-width',dest='wrap_width',type='int',default=50,help='wrap output sequence at this number of bases, 0 indicates no wrap (sequence ends up on single line) [default: %default]')
parser.add_option('--batch',dest='batch',action='store_true',help='run in batch mode, interpret arguments after nib file as queries')
parser.add_option('--batch-format',dest='batch_format',type='choice',choices=['BED'],default='BED',help='format to interpret batch files [default: %default]')
#parser.add_option('--mask-type',dest='mask_type',type='choice',choices=['NOMASK','MASK','HARDMASK'],default='NOMASK',help='how to handle masked positions, correspond to original nibFrag options --masked and --hardMasked [default: %default]')

# original nibFrag usage:
#nibFrag - Extract part of a nib file as .fa (all bases/gaps lower case by default)
#usage:
#   nibFrag [options] file.nib start end strand out.fa
#where strand is + (plus) or m (minus)
#options:
#   -masked - use lower case characters for bases meant to be masked out
#   -hardMasked - use upper case for not masked-out and 'N' characters for masked-out bases
#   -upper - use upper case characters for all bases
#   -name=name Use given name after '>' in output sequence
#   -dbHeader=db Add full database info to the header, with or without -name option
#   -tbaHeader=db Format header for compatibility with tba, takes database name as argument

# original nibFrag options
nibFrag_grp = OptionGroup(parser,"Original nibFrag options")
nibFrag_grp.add_option('--masked',dest='masked',action='store_true',help='use lower case characters for bases meant to be masked out')
nibFrag_grp.add_option('--hardMasked',dest='hardmasked',action='store_true',help='use upper case for non masked-out and \'N\' characters for masked-out bases')
nibFrag_grp.add_option('--upper',dest='upper',action='store_true',help='use upper case characters for all bases')
nibFrag_grp.add_option('--name',dest='name',default=None,help='Use given name after \'>\' in output sequence')
nibFrag_grp.add_option('--dbHeader',dest='dbHeader',default=None,help='Add full database info to the header, with or without -name option')
nibFrag_grp.add_option('--tbaHeader',dest='tbaHeader',default=None,help='Format header for compatibility with tba, takes database name as argument')
parser.add_option_group(nibFrag_grp)


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 1 :
        parser.print_usage()
        parser.exit(1)

    # setup
    nib_path = args[0]
    nib_dir,nib_fn,nib_base,nib_ext = get_file_parts(nib_path)

    queries = []
    if opts.batch :

        if len(args) < 2 :
            parser.error('Two arguments must be supplied in batch mode')

        batch_fns = args[1:]

        for fn in batch_fns :
            if opts.batch_format == 'BED' :
                for bed in BEDFile(fn) :
                    if bed['chrom'] != nib_base :
                        warnings.warn('Chromosome in BED line %s does not match file %s, skipping'%(bed['chrom'],nib_base))
                    else :
                        queries.append((int(bed['chromStart']),int(bed['chromEnd']),bed['strand']))
    else :

        if len(args) < 4 :
            parser.error('Four arguments must be supplied in non-batch mode')

        # setup
        strand = args[3]
        start, end = int(args[1]),args[2]
        if end == 'END' :
            end = -1
        else :
            end = int(end)
            if end < start :
                parser.error('Stop coordinate %d smaller than start %d'%(end,start))

        queries.append((start,end,strand))

    mask_type = NOMASK
    if opts.masked :
        mask_type = MASK
    elif opts.hardmasked :
        mask_type = HARDMASK

    # set the output file
    if len(args) > 4 :
        out_f = open(args[4],'w')
    else :
        out_f = sys.stdout

    # get the sequences from the .nib file
    try :
        seqs = get_nib_seq_batch(nib_path,queries,mask_type)
    except NibException, e :
        sys.stderr.write(e.message+'\n')
        sys.exit(1)

    # construct header
    nbases = validate_nib_file(nib_path)

    # output all queries
    for query, seq in zip(queries,seqs) :
        start,end,strand = query
        if end == -1 :
            end = nbases
        fields = {}
        fields['name'] = nib_path+':%d-%d'%(start,end) if not opts.name else opts.name
        fields['db'] = ''

        if opts.tbaHeader :
            # ignored for some reason in nibFrag when tbaHeader supplied and dbHeader is not
            fields['name'] = '' if not opts.dbHeader else fields['name']
            fields['db'] = '%s.%s:%d-%d of %d'%(opts.tbaHeader,nib_base,start,end,nbases)
        if opts.dbHeader :
            fields['db'] = ':%s.%s:%d-%d:%s:%d'%(opts.dbHeader,nib_base,start,end,strand,nbases)

        header_tmpl = '>%(name)s%(db)s\n'

        # write output
        out_f.write(header_tmpl%fields)
        if opts.upper :
            seq = seq.upper()
        if opts.wrap_width == 0 :
            out_f.write(seq+'\n')
        else :
            for i in xrange(0,len(seq),opts.wrap_width) :
                out_f.write(seq[i:i+opts.wrap_width]+'\n')

