#!/usr/bin/env python

import os
import sys
from csv import DictReader, DictWriter, QUOTE_NONE
from optparse import OptionParser

from chipsequtil import KnownGeneFile, get_file_parts

#args = ['/nfs/genomes/mouse_gp_jul_07/anno/knownGene-2010-07-08.txt','/nfs/genomes/mouse_gp_jul_07/anno/kgXref-2010-07-08.txt']
args = ['/nfs/genomes/mouse_gp_jul_07/anno/knownGene-2010-08-03.gtf','/nfs/genomes/mouse_gp_jul_07/anno/kgXref-2010-07-08.txt']
usage = '%prog <knownGene annotation>'
description = 'convert a UCSC knownGene annotation to GFF'
parser = OptionParser(usage=usage,description=description)


if __name__ == '__main__' :

    opts, args = parser.parse_args(args)

    kg_path,kg_fn,kg_base,kg_ext = get_file_parts(args[0])
    #kg_f = KnownGeneFile(args[0])

    # xref for finding gene symbols
    kgXref_fn = args[1]
    kgXref_fieldnames = ['kgID','mRNA','spID','spDisplayID','geneSymbol','refseq','proAcc','description']
    xref_map = dict([(x['kgID'],x) for x in DictReader(open(kgXref_fn),delimiter='\t',fieldnames=kgXref_fieldnames)])

    gff_headers = ['seqname','source','feature','start','end','score','strand','frame','attributes']
    gff_reader = DictReader(open(args[0]),delimiter='\t',fieldnames=gff_headers)
    gff_writer = DictWriter(sys.stdout,delimiter='\t',fieldnames=gff_headers,quotechar='',quoting=QUOTE_NONE,lineterminator='\n')
    #gff_writer.writerow(dict([(x,x) for x in gff_headers]))

    for i,rec in enumerate(gff_reader) :
        #d = {}
        #d['seqname'] = rec['chrom']
        #d['source'] = 'UCSC_knownGene'
        #d['feature'] = 'gene'
        #d['start'] = rec['txStart']
        #d['end'] = rec['txEnd']
        #d['score'] = '.'
        #d['strand'] = rec['strand']
        #d['frame'] = '.'
        #gene_name = rec['name']

        gff_attrs_lst = [x.strip() for x in rec['attributes'].split(';')][:-1]
        gff_attrs = {}
        for attr in gff_attrs_lst :
            k,v = attr.split(' ',1)
            gff_attrs[k] = eval(v)

        kg_name = gff_attrs['gene_id']

        # try to find a gene symbol
        gene_id = xref_map[kg_name].get('geneSymbol',None)
        #gene_id = kg_name
        #if gene_id is None :
        #    gene_id = xref_map[kg_name].get('mRNA',None)
        #if gene_id is None :
        #    gene_id = xref_map[kg_name].get('refseq',None)
        if gene_id is None : # I give up
            gene_id = kg_name

        gff_attrs_lst += ['gene_name "%s"'%gene_id]
        rec['attributes'] = '; '.join(gff_attrs_lst)
        gff_writer.writerow(rec)

        # now write the exons
        #d['feature'] = 'exon'
        #for j,(st,en) in enumerate(zip(rec['exonStarts'],rec['exonEnds'])) :
        #    d['start'] = st
        #    d['end'] = en
        #    d['attributes'] = '; '.join(['gene_id "%s"'%gene_id,'transcript_id "%s"'%rec['name'],'exon_number "%d"'%(j+1),'ID "%s.exon_%d"'%(rec['name'],j),'PARENT "%s"'%rec['name']])
        #    gff_writer.writerow(d)


    # version with knownGene in gene_name
    # version with symbol in gene_name
