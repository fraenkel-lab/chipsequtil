#!/usr/bin/env python

import sys

from optparse import OptionParser

from chipsequtil import check_org_settings
from chipsequtil.util import MultiLineHelpFormatter
from chipsequtil.sampling import rejection_sample_bg
from TAMO.seq import Fasta

usage = '%prog [options] <organism> <fasta file> [<fasta file> ... ]'
description = """Use rejection sampling to generate a set of background/random
sequences matching the distance to nearest transcription start site, sequence
length, and GC content distributions of the input fasta file(s).  Generated
sequences are genomic sequences sampled based on these distributions. All sequences
from all files are used to generate the background sequences. The following
command must output a path to a nib genomic sequence directory and refGene
annotation, respectively :

$> org_settings.py <organism> genome_dir
$> org_settings.py <organism> refgene_anno_path

Utility prints out generated fasta records to stdout by default.
"""
epilog = "Note: script only considers sequences with unique header names, only the last record of those with identical header names is used"
parser = OptionParser(usage=usage,description=description,formatter=MultiLineHelpFormatter())
parser.add_option('-n','--num-seqs',dest='num_seqs',default='1x', help='number of sequences to generate, either absolute number or factor of # input sequences, e.g. 2.5x for 2.5 times the # of input sequences [default: 1x]')
parser.add_option('--output',dest='output',default=None,help='file to output fasta records to [default: stdout]')
parser.add_option('--bed',dest='bed',action='store_true', help='also produce a BED formatted file representing sampled sequences')
parser.add_option('--bed-output',dest='bed_output',default='output.bed',help='with --bed, file to output BED records to [default: %default]')

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 2 :
        parser.error('Must be 2 non-option arguments')

    organism, fasta_fns = args[0], args[1:]

    reqd_settings = ['genome_dir','refgene_anno_path']
    if not check_org_settings(organism,reqd_settings) :
        parser.error('The <organism> settings set must contain paths for %s'%reqd_settings)

    # load up all the fasta records
    fasta_recs = {}
    for fasta_fn in fasta_fns :
        fasta = Fasta.load(fasta_fn)
        fasta_recs.update(fasta)

    # parse --num-seqs argument
    if opts.num_seqs.endswith('x') :
        num_seq_factor = float(opts.num_seqs[:-1])
        num_seqs = int(len(fasta_recs)*num_seq_factor)
    else :
        try :
            num_seqs = int(opts.num_seqs)
        except TypeError :
            parser.error("Incorrect format of --num-seqs argument, must either be an integer or a factor ending with x, e.g. 2.5x")

    # generate the sequences
    gen_seqs = rejection_sample_bg(fasta_recs,organism,num_samples=num_seqs)

    # write out to file
    if opts.output :
        Fasta.write(gen_seqs,opts.output)
    else :
        sys.stdout.write(''.join(['>%s\n%s\n'%(k,v) for k,v in gen_seqs.items()]))

    if opts.bed :
        bed_f = open(opts.bed_output,'w')
        bed_f.write(''.join([k.replace(':','\t').replace('-','\t')+'\n' for k in gen_seqs.keys()]))
        bed_f.close()

