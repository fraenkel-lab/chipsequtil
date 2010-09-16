#!/usr/bin/env python

import gzip
import sys
from collections import defaultdict as dd
from csv import DictReader, DictWriter
from optparse import OptionParser
from sqlite3 import connect

from chipsequtil import KnownGeneFile

# TODO make these parameters?
#affy_anno_fn = 'Mouse430A_2.na30.annot.csv'

usage = '%prog [options] <knownGene annotation> <knownToMOE430 file> <knownGene Xref file> <microarray data file>'
description = 'Maps probset data to knownGene database provided by UCSC. Probesets \
that map to multiple knownGenes have one record per knownGene with duplicate data \
otherwise.  Output is knownGene id prepended to each record in microarray data file.'
parser = OptionParser(usage=usage,description=description)
parser.add_option('--output',dest='output',default=None,help='file to output mapping to [default: stdout]')
#parser.add_option('--symbol-xref',dest='symbol_xref',default=None,help='use the provided kgXref file to output gene symbols as second column')

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    #affy_bioc_fn = 'microarray_analysis/cbfb_vector_BH_all.txt'
    #knownToMOE_sql_fn = 'knownToMOE430.sql'
    #knownToMOE_data_fn = 'knownToMOE430.txt'

    if len(args) < 3 :
        parser.error('Incorrect number of arguments provided')

    known_gene_fn = args[0]
    knownToMOE_data_fn = args[1]
    Xref_fn = args[2]
    affy_bioc_fn = args[3]

    # affymetrix file from bioconductor
    affy_bioc_f = open(affy_bioc_fn)
    affy_bioc = {}
    affy_bioc_reader = DictReader(affy_bioc_f,delimiter="\t")
    for row in  affy_bioc_reader :
        affy_bioc[row['ID']] = row

    # knownGene annotation
    kg = KnownGeneFile(known_gene_fn)
    kg_ids = dict([(x['name'],x) for x in kg])

    # affy to knownGene
    affy_to_kg_map = dd(list)
    affy_to_kg_fields = ['kgID','affyID']
    affy_to_kg_f = open(knownToMOE_data_fn)
    kg_to_affy_map = dd(list)
    for row in DictReader(affy_to_kg_f,fieldnames=affy_to_kg_fields,delimiter="\t") :
        affy_to_kg_map[row['affyID'][2:]].append(row['kgID'])
        kg_to_affy_map[row['kgID']].append(row['affyID'][2:])

    if opts.output :
        out_f = open(opts.output,'w')
    else :
        out_f = sys.stdout

    out_header = ['knownGeneID']+affy_bioc_reader.fieldnames

    # see if the user wants gene symbols too
    opts.symbol_xref = Xref_fn
    if opts.symbol_xref :
        kgXref_fieldnames = ['kgID','mRNA','spID','spDisplayID','geneSymbol','refseq','protAcc','description']
        symbol_xref_reader = DictReader(open(opts.symbol_xref),fieldnames=kgXref_fieldnames,delimiter='\t')
        symbol_xref_map = {}
        for rec in symbol_xref_reader :
            symbol_xref_map[rec['kgID']] = rec
        out_header = ['knownGeneID','geneSymbol']+affy_bioc_reader.fieldnames

    out_writer = DictWriter(out_f,delimiter='\t',fieldnames=out_header,lineterminator='\n')
    out_writer.writerow(dict(zip(out_header,out_header)))
    for probesetID, data in affy_bioc.items() :
        kg_ids = affy_to_kg_map[probesetID]
        for kg_id in kg_ids :
            out_l = {'knownGeneID':kg_id}
            if opts.symbol_xref :
                out_l['geneSymbol'] = symbol_xref_map[kg_id]['geneSymbol']
            out_l.update(data)
            out_writer.writerow(out_l)

    # figure out if any probsets map to non-overlapping loci
    # dirty dirty dirty dirty
    if False :
        affy_id_loci = {}
        for affyID, kgIDs in affy_to_kg_map.items() :
            # check all pairwise kgIDs to make sure they all overlap in transcription start sites
            kg_id_loci = dd(list)
            for i, kgID1 in enumerate(kgIDs) :
                kgID1_rec = kg_ids[kgID1]
                kg_id_loci[kgID1].append(kgID1_rec)
                for j, kgID2 in enumerate(kgIDs) :
                    kgID2_rec = kg_ids[kgID2]
                    # these are all gene overlap conditions
                    #kg1Start = kgID1_rec['txEnd'] if kgID1_rec['strand'] == '-' else kgID1_rec['txStart']
                    #kg1End = kgID1_rec['txStart'] if kgID1_rec['strand'] == '-' else kgID1_rec['txEnd']
                    #kg2Start = kgID2_rec['txEnd'] if kgID2_rec['strand'] == '-' else kgID2_rec['txStart']
                    #kg2End = kgID2_rec['txStart'] if kgID2_rec['strand'] == '-' else kgID2_rec['txEnd']
                    kg1Start, kg1End = kgID1_rec['txStart'], kgID1_rec['txEnd']
                    kg2Start, kg2End = kgID2_rec['txStart'], kgID2_rec['txEnd']
                    if (kg2Start <= kg1Start <= kg2End or \
                       kg1Start <= kg2Start <= kg1End or \
                       (kg2Start < kg1Start and kg2End > kg1End) or \
                       (kg1Start < kg2Start and kg1End > kg2End)) and \
                       kgID1_rec['chrom'] == kgID2_rec['chrom'] and \
                       i != j :
                        # we have overlap
                        pass
                    elif i != j :
                        # doesn't overlap oh noes
                        kg_id_loci[kgID1].append(kgID2_rec)
            for kg_id, kg_recs in kg_id_loci.items() :
                if len(kg_recs) != 1 :
                    affy_id_loci[affyID] = (kg_id, len(kg_recs),len(kgIDs),kg_recs,kgIDs)

        if len(affy_id_loci) != 0 :
            sys.stderr.write('Warning: %d probeset ids map to non-overlapping loci'%len(affy_id_loci))


