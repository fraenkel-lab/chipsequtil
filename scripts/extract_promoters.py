#!/usr/bin/env python

import sys
from csv import writer
from optparse import OptionParser


from chipsequtil import get_org_settings, RefGeneFile
from chipsequtil.nib import NibDB

parser = OptionParser()
parser.add_option('-u','--upstream',type='int',default=3000,help='upstream window from TSS to extract')
parser.add_option('-d','--downstream',type='int',default=1000,help='downstream window from TSS to extract')

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    org_settings = get_org_settings(args[0])

    refgene_fn = org_settings['refgene_anno_path']
    refgene_f = RefGeneFile(refgene_fn)

    nib_db = NibDB(nib_dirs=[org_settings['genome_dir']])

    out_f = writer(sys.stdout,delimiter='\t')
    seq_recs = []
    for rec in refgene_f :
        st, end = max(0,int(rec['txStart'])-opts.upstream), min(int(rec['txStart'])+opts.downstream,nib_db.db_info[rec['chrom']]['nbases'])
        seq_recs.append((rec['chrom'],st,end,rec['strand']))

    fasta_recs = nib_db.get_fasta_batch(seq_recs)

    for header, seq in zip(*fasta_recs) :
        sys.stdout.write(header+seq+'\n')
