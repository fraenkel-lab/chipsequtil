#!/usr/bin/env python

import sys
import warnings

from csv import reader, writer
from collections import defaultdict as dd
from optparse import OptionParser

usage = '%prog -b <mapped DNA binding file>|-a <mapped microarray file> [-b <mapped DNA binding file> ...] [-a <mapped microarray file> ...]'
description = """Join all files on the first column, concatenating records with \
matching entries onto one line per entry.  Understands DNA binding data as mapped \
with *map_peaks_to_known_genes.py* utility microarray data as mapped by \
*probeset_to_known_genes.py* utility, passed to program using *-b* and *-a* options \
respectively.  If a file contains more than one mapping to a gene additional columns \
are added. At least one file of either type is required.  Field names are written as \
<filename>.<original field name>.<map number>
"""
epilog="Note: microarray files should have been created by bioconductor, and all files should have a row of fieldnames as the first line"
parser = OptionParser(usage=usage,description=description,epilog=epilog)
parser.add_option('-a','--affy-file',dest='affy_file',action='append',default=[],help='add a mapped microarray file')
parser.add_option('-b','--bind-file',dest='bind_file',action='append',default=[],help='add a mapped DNA binding file (e.g. MACS, BED)')
#parser.add_option('-b','--bed-file',dest='bed_file',action='append',default=[],help='add a mapped BED formatted peaks file')
parser.add_option('-m','--macs-file',dest='macs_file',action='append',default=[],help='DEPRECATED: use -b instead, add a mapped default MACS formatted peaks (*.xls) file')
parser.add_option('--output',dest='output',default=None,help='file to output joined records to [default: stdout]')
#parser.add_option('--intersect',dest='intersect',action='store_true',help='only output records common to all file passed in')
parser.add_option('--first-only',dest='first_only',action='store_true',help='only output the first mapping to a gene from each file')
parser.add_option('--binary',dest='binary',action='store_true',help='output only one column per file with a 0 or 1 to indicate whether a mapping exists in that file')
parser.add_option('--binary-plus',dest='binary_plus',action='store_true',help='output one column per file with a 0 or 1 to indicate whether a mapping exists in that file in addition to all other columns')
parser.add_option('--field-types',dest='field_types',action='store_true',help='prepend BIND or AFFY to the beginning of all appropriate columns')
#parser.add_option('--symbols',dest='symbols',action='store_true',help='mapped files contain symbols in second column (per map_peaks_to_known_genes.py|probeset_to_known_gene.py --symbol-xref option)')

if __name__ == '__main__' :

    opts,args = parser.parse_args(sys.argv[1:])

    if len(args) > 0 :
        parser.error('There were non-option command line arguments passed, all files should have a preceeding option indicating filetype')

    if len(opts.macs_file) != 0 :
        warnings.warn('The -m option is deprecated, please replace these flags with -b instead.  Adding MACS filenames to binding filename list.',DeprecationWarning)
        opts.bind_file.extend(opts.macs_file)

    if len(opts.bind_file) == 0 and len(opts.affy_file) == 0 :
        parser.error('No files were passed in, aborting')

    # union of all genes
    all_genes = set()

    # TODO - fix intersect w/ binary
    opts.intersect = False

    # TODO - actually make this an option, or the default
    opts.symbols = True
    if opts.symbols :
        symbol_map = {}

    # read all the files in
    def get_file_dict(fns,header_prefix='') :
        file_map = dd(lambda: dd(list))
        out_fieldnames = []
        blank_entry = []
        for fn in fns :
            max_maps = 0
            f = reader(open(fn),delimiter='\t')
            #f = open(fn)
            fieldnames = f.next()
            fieldnames = fieldnames[2:] # we don't want existing knownGeneID or geneSymbol
            # read in the data, create a dictionary
            for l in f :
                if opts.symbols :
                    gene, symbol, data = l[0],l[1],l[2:]
                    symbol_map[gene] = symbol
                else :
                    gene, data = l.split('\t',1)
                file_map[fn][gene].append(data)
                max_maps = max(max_maps,len(file_map[fn][gene]))
                all_genes.add(gene)

            # if we're adding a binary column, do it
            if opts.binary_plus :
                out_fieldnames.append(header_prefix+fn+'.MAPPED')

            # construct the fieldnames for this file
            for i in range(max_maps) :
                out_fieldnames.extend(['%s%s.%d.%s'%(header_prefix,fn,i,h) for h in fieldnames])

            # pad out data entries w/ fewer than max_maps
            for gene,data in file_map[fn].items() :
                while len(data) < max_maps :
                    data.append(['']*len(fieldnames))
            file_map[fn]['blank'] = [['']*len(fieldnames) for _ in range(max_maps)]
        return file_map,out_fieldnames

    #macs_file_map, macs_fieldnames = get_file_dict(opts.macs_file)
    #bed_file_map, bed_fieldnames = get_file_dict(opts.bed_file)
    bind_prefix = 'BIND:' if opts.field_types else ''
    affy_prefix = 'AFFY:' if opts.field_types else ''
    bind_file_map, bind_fieldnames = get_file_dict(opts.bind_file,bind_prefix)
    affy_file_map, affy_fieldnames = get_file_dict(opts.affy_file,affy_prefix)

    # prepare output objects
    out_f = open(opts.output,'w') if opts.output else sys.stdout
    map_fieldnames = ['knownGeneID']
    if opts.symbols :
        map_fieldnames.append('geneSymbol')
    #all_fieldnames = map_fieldnames+macs_fieldnames+bed_fieldnames+affy_fieldnames
    all_fieldnames = map_fieldnames+bind_fieldnames+affy_fieldnames
    if opts.binary :
        #all_fieldnames = map_fieldnames+opts.macs_file+opts.bed_file+opts.affy_file
        all_fieldnames = [x+'.MAPPED' for x in map_fieldnames+opts.bind_file+opts.affy_file]
    join_writer = writer(out_f,delimiter='\t')
    join_writer.writerow(all_fieldnames)

    # go through all the genes and print out lines
    for gene in all_genes :
        gene_line = [gene]
        if opts.symbols :
            gene_line.append(symbol_map[gene])
        #for filetype_data,fns in zip([macs_file_map,bed_file_map,affy_file_map],[opts.macs_file,opts.bed_file,opts.affy_file]) :
        for filetype_data,fns in zip([bind_file_map,affy_file_map],[opts.bind_file,opts.affy_file]) :
            for fn,recs in [(fn,filetype_data[fn]) for fn in fns] :
            #for fn,recs in d.items() :
                if recs.has_key(gene) :
                    # only output the first entry
                    if opts.first_only :
                        gene_line.extend(recs[gene][0])
                    # only output a 1 or a zero
                    elif opts.binary :
                        gene_line.extend('1')
                    # else output normally
                    else :
                        # add binary column in addition to other output
                        if opts.binary_plus :
                            gene_line.extend('1')
                        for rec in recs[gene] :
                            gene_line.extend(rec)
                else :
                    # if intersecting, ignore this gene
                    if opts.intersect :
                        continue
                    elif opts.binary :
                        gene_line.extend('0')
                    else :
                        # add binary column in addition to other output
                        if opts.binary_plus :
                            gene_line.extend('0')
                        for blank in filetype_data[fn]['blank'] :
                            #print len(blank)
                            gene_line.extend(blank)
                #print fn, gene_line[2], len(gene_line), gene_line
        join_writer.writerow(gene_line)

    if opts.output : out_f.close()
