#!/usr/bin/env python

import re
import sys
from csv import writer
from optparse import OptionParser

from collections import defaultdict

from chipsequtil import get_org_settings, RefGeneFile
from chipsequtil.nib import NibDB
from chipsequtil.util import MultiLineHelpFormatter as MF

usage = "%prog [options] <organism>"
description = """Extract the promoter sequences in FASTA format from all genes
or a list of genes specified in an input file.  Gene annotation is RefGene
corresponding to the organism passed in, paths returned by:

$> org_settings.py <organism> refgene_anno_path
$> org_settings.py <organism> genome_dir

must be valid."""
parser = OptionParser(usage=usage,description=description,formatter=MF())
parser.add_option('-u','--upstream',type='int',default=3000,help='upstream window from TSS to extract [default: %default]')
parser.add_option('-d','--downstream',type='int',default=1000,help='downstream window from TSS to extract [default: %default]')
parser.add_option('-l','--gene-list',dest='gene_list',default=None,
                  help='file containing a list of gene identifiers to extract, one per line [default: %default]')
gene_type_choices = ['symbol','refgene']
parser.add_option('-t','--gene-type',dest='gene_type',type='choice',
                  choices=gene_type_choices,default=gene_type_choices[0],
                  help='type of gene identifier in gene list, choose from %s [default: %%default]'%gene_type_choices)
parser.add_option('-o','--output',dest='output',default=None,
                  help='file to write fasta records to [default: stdout]')

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 1 :
        parser.error('Exactly one argument is required')

    org_settings = get_org_settings(args[0])

    refgene_fn = org_settings['refgene_anno_path']
    refgene_f = RefGeneFile(refgene_fn)

    nib_db = NibDB(nib_dirs=[org_settings['genome_dir']])

    gene_list = None
    if opts.gene_list :
        gene_list = [x.strip() for x in open(opts.gene_list).readlines()]

    id_index = 'bin'
    if opts.gene_type != gene_type_choices[0] :
        if opts.gene_type  == 'refgene' :
            id_index = 'name'

    seq_recs = []
    gene_map = defaultdict(list)
    for rec in refgene_f :
        if gene_list and rec[id_index] not in gene_list : continue # skip this one
        st, end = max(0,int(rec['txStart'])-opts.upstream), min(int(rec['txStart'])+opts.downstream,nib_db.db_info[rec['chrom']]['nbases'])
        key = (rec['chrom'],st,end,rec['strand'])
        seq_recs.append(key)
        gene_map[key[:-1]].append(rec['bin']+'/'+rec['name'])

    fasta_recs = nib_db.get_fasta_batch(seq_recs)

    out_f = open(opts.output,'w') if opts.output else sys.stdout
    header_regex = re.compile('^.*(chr[0-9MXY]+).*:([0-9]+)-([0-9]+).*$')
    for header, seq in zip(*fasta_recs) :
        # map sequences back to gene names using the header
        reg_obj = header_regex.search(header)
        if reg_obj is not None :
            chrm,st,end = reg_obj.groups()
            gene_names = gene_map.get((chrm,int(st),int(end)))
            if gene_names is not None :
                header = header.strip()+':'+','.join(gene_names)+'\n'
        out_f.write(header+seq+'\n')
