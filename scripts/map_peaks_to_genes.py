#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from collections import defaultdict as dd
from chipsequtil import MACSOutput, BEDOutput, RefGeneOutput, parse_number
from csv import DictReader, DictWriter

usage = '%prog [options] <refGene file> <peaks file>'
description = """
Map the peaks in <peaks file> to genes in <refGene file>.  <refGene file> is
format is as specified in http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.sql.
<peaks file> format is as produced by MACS."""
epilog = ''
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('--upstream-window',dest='upst_win',type='int',default=5500,help='window width in base pairs to consider promoter region [default: %default]')
parser.add_option('--downstream-window',dest='dnst_win',type='int',default=2500,help='window width in base pairs to consider downstream region [default: %default]')
parser.add_option('--map-output',dest='peak_output',default=sys.stdout,help='filename to output mapped peaks in BED format to [default: stdout]')
parser.add_option('--stats-output',dest='stats_output',default=sys.stderr,help='filename to output summary stats in conversion [default: stderr]')
parser.add_option('--peaks-format',dest='peaks_fmt',default='MACS',type='choice',choices=['MACS','BED'],help='format of peaks input file [default: %default]')

# TODO - options
#parser.add_option('--use-cds',dest='use_cds',action='store_true',help='use cdsStart and cdsEnd fields instead of txStart and txEnd to do mapping')
#parser.add_option('--capture-intergenic'...)
#parser.add_option('--map-format',dest='peak_format',type='choice',choices=['default','BED'],help='format of peak output [default: %default]')
#parser.add_option('--stats-format',dest='stats_format',type='choice',choices=['human','python'],help='format of summary stats output [default: %default]')

def parse_gene_ref(ref_gene) :
    #FIXME - maybe, if galaxy doesn't work out, figure out how to deal with multiple RefGene mapping formats?
    fieldnames = ['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
    reader = DictReader(ref_gene,fieldnames=fieldnames,delimiter='\t')
    gene_ref = dd(list)
    for ref_dict in reader :
        for k,v in ref_dict.items() :
            # coerce numbers where possible
            ref_dict[k] = parse_number(v)

        # turn 'x,x,x,...' into a list
        ref_dict['exonStarts'] = [parse_number(x) for x in ref_dict['exonStarts'].split(',')]
        if ref_dict['exonStarts'][-1] == '' : ref_dict['exonStarts'].remove('')
        ref_dict['exonEnds'] = [parse_number(x) for x in ref_dict['exonEnds'].split(',')]
        if ref_dict['exonEnds'][-1] == '' : ref_dict['exonEnds'].remove('')

        gene_ref[ref_dict['chrom']].append(ref_dict)

    return gene_ref


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 2 :
        parser.error('Must provide two filename arguments')

    gene_ref = parse_gene_ref(open(args[0]))
    if opts.peaks_fmt == 'MACS' :
        fieldnames = MACSOutput.FIELD_NAMES
        chr_field, start_field, end_field = 'chr', 'start', 'end'
    elif opts.peaks_fmt == 'BED' :
        fieldnames = BEDOutput.FIELD_NAMES
        chr_field, start_field, end_field = 'chrom', 'chromStart', 'chromEnd'
    else :
        fieldnames = []

    peaks_reader = DictReader(open(args[1]),fieldnames=fieldnames,delimiter='\t')

    # default output format:
    # <chromo> <peak loc> <accession #> <gene symbol> <strand> <map type> <map subtype> <score> <dist from feature>
    # score = (peak-TSS)/(TSE-TSS) - peak distance from TSS as fraction of length of gene
    output_fields = ['chromo',
                     'peak loc',
                     'accession #',
                     'gene symbol',
                     'strand',
                     'map type',
                     'map subtype',
                     'score',
                     'dist from feature',
    ]
    if opts.peak_output != sys.stdout :
        opts.peak_output = open(opts.peak_output,'w')
    peaks_writer = DictWriter(opts.peak_output,output_fields,delimiter='\t',lineterminator='\n')
    unique_genes = set()
    map_stats = dd(int)
    for peak in peaks_reader :

        # if this is a comment or header line get skip it
        if peak[fieldnames[0]].startswith('#') or \
           peak[fieldnames[0]] == fieldnames[0] or \
           peak[fieldnames[0]].startswith('track') : continue

        # coerce values to numeric if possible
        for k,v in peak.items() : peak[k] = parse_number(v)

        # peak assumed to be in the middle of the reported peak range
        peak_loc = (peak[start_field]+peak[end_field])/2

        chrom_genes = gene_ref[peak[chr_field]]

        if len(chrom_genes) == 0 :
            sys.stderr.write('WARNING: peak chromosome %s not found in gene reference, skipping: %s\n'%(peak[chr_field],peak))
            continue

        mapped = False

        # walk through the genes for this chromosome
        for gene in chrom_genes :

            # reusable dictionary for output
            out_d = {}.fromkeys(output_fields,0)
            out_d['map type'] = ''
            out_d['chromo'] = peak[chr_field]
            out_d['peak loc'] = peak_loc

            # determine intervals for promoter, gene, and downstream
            if gene['strand'] == '+' :
                promoter_coords = max(gene['txStart']-1-opts.upst_win,0), gene['txStart']-1
                gene_coords = gene['txStart'], gene['txEnd']
                downstream_coords = gene['txEnd']+1, gene['txEnd']+1+opts.dnst_win
            else :
                promoter_coords = gene['txEnd']+1, gene['txEnd']+1+opts.upst_win # +1 because we're using 1 based indexing
                gene_coords = gene['txStart'], gene['txEnd']
                downstream_coords = gene['txStart']-1-opts.dnst_win, gene['txStart']-1 # -1 because we're using 1 based indexing

            # check for promoter
            if peak_loc >= promoter_coords[0] and peak_loc <= promoter_coords[1] :
                out_d['map type'] = 'promoter'
                out_d['dist from feature'] = peak_loc - promoter_coords[1] if gene['strand'] == '+' else promoter_coords[0] - peak_loc

            # check for gene
            elif peak_loc >= gene_coords[0] and peak_loc <= gene_coords[1] :
                # check for intron/exon
                exon_coords = zip(gene['exonStarts'],gene['exonEnds'])
                in_exon = False
                for st,en in exon_coords :
                    if peak_loc >= st and peak_loc <= en :
                        in_exon = True
                        break
                out_d['map type'] = 'gene'
                out_d['map subtype'] = 'exon' if in_exon else 'intron'

                # score = (peak-TSS)/(TSE-TSS) - peak distance from TSS as fraction of length of gene
                gene_len = float(gene_coords[1]-gene_coords[0])
                out_d['score'] = (peak_loc-gene_coords[0])/gene_len if gene['strand'] == '+' else (gene_coords[1]-peak_loc)/gene_len

                # distance calculated from start of gene
                out_d['dist from feature'] = peak_loc - promoter_coords[1] if gene['strand'] == '+' else promoter_coords[0] - peak_loc

                map_stats[out_d['map subtype']] += 1

            # check for downstream
            elif peak_loc >= downstream_coords[0] and peak_loc <= downstream_coords[1] :
                out_d['map type'] = 'after'
                out_d['dist from feature'] = peak_loc - downstream_coords[0] if gene['strand'] == '+' else downstream_coords[1] - peak_loc

            # does not map to this gene
            else :
                pass

            # map type is not blank if we mapped to something
            if out_d['map type'] != '' :

                out_d['accession #'] = gene['name']
                out_d['gene symbol'] = gene['geneName']
                out_d['strand'] = gene['strand']

                map_stats[out_d['map type']] += 1
                peaks_writer.writerow(out_d)

                unique_genes.add(gene['name'])
                mapped = True

                """
                print 'Peak:',peak
                print 'Gene:',gene
                print 'Peak loc:',peak_loc
                print promoter_coords
                print gene_coords
                print downstream_coords
                raw_input('Wait for it...')
                """

                # reset map_type
                out_d['map type'] = ''

        if not mapped :
            #out_d['map type'] = 'intergenic'
            #peaks_writer.writerow(out_d)
            map_stats['intergenic'] += 1

    if opts.peak_output != sys.stdout :
        opts.peak_output.close()

    if opts.stats_output != sys.stderr :
        opts.stats_output = open(opts.stats_output,'w')

    for k,v in map_stats.items() :
        opts.stats_output.write('%s: %s\n'%(k,v))

    if opts.stats_output != sys.stderr :
        opts.stats_output.close()
