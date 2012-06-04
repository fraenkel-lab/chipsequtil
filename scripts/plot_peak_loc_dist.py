#!/usr/bin/env python

import matplotlib
matplotlib.use('AGG')

import matplotlib.pyplot as mp
import numpy as np
import os
import sys

from collections import defaultdict
from csv import reader, writer
from optparse import OptionParser
from StringIO import StringIO

from chipsequtil import MACSFile, BEDFile


usage = '%prog [options] <peaks fn> <gene list fn>'
desc = """Produce a pie chart of the locations of peaks in different bins
(promoter, gene, exon, intron, etc.) and, optionally, save the different
records to their own files for subsequent analysis.  Also produce a histogram
of distance from feature values in mapping file. Peaks file is expected
to be as output by MACS, or alternately as a BED file but then the -b plot
is not available.  Gene list file is expected to be in the format as
output by peaks_to_known_genes.py script."""
parser = OptionParser(usage=usage,description=desc)
parser.add_option('-b','--bar-fn',dest='bar_fn',default=None,help='filename for pvalue stacked bar chart')
parser.add_option('-g','--gene-pie-fn',dest='gene_pie_fn',default=None,help='filename for pie chart image')
parser.add_option('-p','--peak-pie-fn',dest='peak_pie_fn',default=None,help='filename for pie chart image')
parser.add_option('-f','--dist-fn',dest='dist_fn',default=None,help='filename for distance from feature image')
parser.add_option('-s','--save',dest='save',action='store_true',help='write out files containing peaks for each category')
parser.add_option('-d','--output-dir',dest='out_dir',default='.',help='output files created by --save option to this directory')
parser.add_option('--no-plot',dest='no_plot',action='store_true',help='dont show (but save) the figure produced')
parser.add_option('--peaks-format',dest='peak_fmt',type='choice',choices=['MACS','BED'],default='MACS',help='format of peaks file, either MACS or BED [default: MACS]')

GENE_FIELD_NAMES = ['knowngene_id','gene_symbol']
LOC_FIELD_NAMES = ['peak_loc','dist_from_feature','score','map_type','map_subtype']
int_or_none = lambda x: int(x) if x != '' else None
float_or_none = lambda x: float(x) if x != '' else None
LOC_FIELD_TYPES = [int_or_none,float_or_none,float_or_none,str,str]


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 2 :
        parser.error('Exactly 2 non-option argument is required')

    peaks_fn, gene_fn = args

    if opts.peak_fmt == 'BED' :
        peaks_f = BEDFile(peaks_fn)
    else :
        peaks_f = MACSFile(peaks_fn)

    gene_reader = reader(open(gene_fn),delimiter='\t')
    gene_recs, macs_recs, loc_recs = [], [], []
    gene_reader.next() # get rid of header

    gene_field_cnt = len(GENE_FIELD_NAMES)
    macs_field_cnt = len(MACSFile.FIELD_NAMES)
    loc_field_cnt = len(LOC_FIELD_NAMES)
    for rec in gene_reader :

        gene_recs.append(dict(zip(GENE_FIELD_NAMES,rec[:gene_field_cnt])))

        # this automatically coerces recs into correct format
        macs_line = []
        for f,x in zip(MACSFile.FIELD_TYPES,rec[gene_field_cnt:gene_field_cnt+macs_field_cnt]) :
            if f is float :
                f = float_or_none
            macs_line.append(f(x))
        macs_recs.append(dict(zip(MACSFile.FIELD_NAMES,macs_line)))

        loc_line = [f(x) for f,x in zip(LOC_FIELD_TYPES,rec[gene_field_cnt+macs_field_cnt:])]
        loc_recs.append(dict(zip(LOC_FIELD_NAMES,loc_line)))

    loc_dist = defaultdict(int)
    unique_peaks = defaultdict(set)
    exon_scores, intron_scores = [], []
    dist_to_features = defaultdict(list)
    pvals = defaultdict(list)

    fn_base, fn_ext = os.path.splitext(gene_fn)
    if opts.save :
        def get_writer(fn) :
            fd = writer(open(fn,'w'),delimiter='\t')
            header = MACSFile.FIELD_NAMES
            if opts.peak_fmt == 'BED' :
                header = BEDFile.FIELD_NAMES
            fd.writerow(GENE_FIELD_NAMES+header+LOC_FIELD_NAMES)
            return fd
        fds = {}

    for gene, peak, loc in zip(gene_recs, macs_recs, loc_recs) :
        # weird case, not sure why this happens
        if loc['map_subtype'] == '0' :
            loc['map_subtype'] = ''
        key = loc['map_type']+'_%s'%loc['map_subtype'] if loc['map_subtype'] != '' else loc['map_type']
        loc_dist[key] += 1
        dist_to_features[key].append(int(loc['dist_from_feature']))
        if opts.peak_fmt == 'MACS' :
            pvals[key].append(float(peak['-10*log10(pvalue)']))

        map_key = '%s:%d-%d'%(peak['chr'],peak['start'],peak['end'])
        unique_peaks[key].add(map_key)

        if key == 'gene_exon' :
            exon_scores.append(loc['score'])
        elif key == 'gene_intron' :
            intron_scores.append(loc['score'])

        if opts.save :
            row = [gene[f] for f in GENE_FIELD_NAMES] + \
                  [peak[f] for f in MACSFile.FIELD_NAMES] + \
                  [loc[f] for f in LOC_FIELD_NAMES]
            if not fds.has_key(key) :
                fn = os.path.join(opts.out_dir,fn_base+'_'+key+fn_ext)
                fds[key] = get_writer(fn)
            fds[key].writerow(row)

    # now find which peaks are intergenic
    intergenic = []
    num_peaks = 0
    all_unique_peaks = reduce(lambda x,y: x.union(y), unique_peaks.values())
    for l in peaks_f :
        peak = l
        map_key = '%s:%d-%d'%(peak['chr'],peak['start'],peak['end'])
        if map_key not in all_unique_peaks :
            unique_peaks['intergenic'].add(map_key)
            intergenic.append(peak)
            if opts.peak_fmt == 'MACS' :
                pvals['intergenic'].append(peak['-10*log10(pvalue)'])
        num_peaks += 1

    num_int = len(intergenic)
    loc_dist['intergenic'] = num_int
    if opts.save :
        fn = os.path.join(opts.out_dir,fn_base+'_intergenic.xls')
        fd = writer(open(fn,'w'),delimiter='\t')
        fd.writerow(MACSFile.FIELD_NAMES)
        fd.writerows([[x[f] for f in MACSFile.FIELD_NAMES] for x in intergenic])

    exon_scores, intron_scores = np.array(exon_scores), np.array(intron_scores)

    font = {'size':'9'}
    mp.rc('font',**font)
    fig = mp.figure(figsize=(4,4))

    bin_order = ('intergenic','gene_exon','gene_intron','promoter','after')
    colors = 'bgrcm'

    # pie chart
    #pie_ax_rect = [0.1,0.35, 0.4125,  0.525 ] # left, bottom, width, height
    pie_ax = fig.add_axes((0.15,0.15,0.7,0.7))
    pie_ax.set_title('Gene map distribution\n%d peaks'%num_peaks)
    pie_labels, pie_values = [], []
    for k in bin_order :
        pie_labels.append(k+'\n%d'%(len(unique_peaks[k])))
        pie_values.append(len(unique_peaks[k]))
    pie_ax.pie(pie_values,labels=pie_labels)

    img_fn = fn_base+'_gene_loc.png' if opts.gene_pie_fn is None else opts.gene_pie_fn
    mp.savefig(img_fn)
    mp.clf()


    fig = mp.figure(figsize=(4,4))
    pie_ax = fig.add_axes((0.15,0.15,0.7,0.7))
    pie_ax.set_title('Peak map distribution\n%d peaks'%num_peaks)
    pie_labels, pie_values = [], []
    for k in bin_order :
        pie_labels.append(k+'\n%d'%(loc_dist[k]))
        pie_values.append(loc_dist[k])
    pie_ax.pie(pie_values,labels=pie_labels)

    img_fn = fn_base+'_peak_loc.png' if opts.peak_pie_fn is None else opts.peak_pie_fn
    mp.savefig(img_fn)
    mp.clf()

    fig = mp.figure(figsize=(4,4))
    # dist to feature histogram
    #hist_ax_rect = [0.65,0.45,0.25,0.45]
    hist_ax = fig.add_axes((0.15,0.15,0.7,0.7))
    hist_ax.set_title('Peak distance from TSS')
    # join all the lists together
    dists = sum(dist_to_features.values(),[])
    pdf, bins, patches = hist_ax.hist(dists,bins=20)
    #h = mp.hist(dists,bins=20)
    hist_ax.set_xlim((int(min(dists)),int(max(dists))))

    dist_fn = fn_base+'_dist.png' if opts.dist_fn is None else opts.dist_fn
    mp.savefig(dist_fn)
    mp.clf()

    if opts.peak_fmt == 'MACS' :
        fig = mp.figure(figsize=(4,4))
        bar_ax = fig.add_axes((0.15,0.15,0.7,0.7))
        pval_hists = {}
        min_pval, max_pval = min([min(v) for v in pvals.values()]), max([max(v) for v in pvals.values()])
        for key,pvals in pvals.items() :
            vals, bins = np.histogram(pvals,range=(0,max_pval),bins=20)
            lv = np.log10(vals)
            lv[np.isneginf(lv)] = 0.1
            pval_hists[key] = lv

        pval_items = [(k,pval_hists[k]) for k in bin_order if pval_hists.has_key(k)]
        bar_width = 0.85*(max_pval-min_pval)/(len(bins)-1)
        print max_pval, min_pval, len(bins)
        print 'bar_width:',bar_width
        bars = []
        b = bar_ax.bar(bins[:-1],pval_items[0][1],width=bar_width,color=colors[0])
        bars.append(b)

        sum_bottoms = pval_items[0][1]
        for i, (key, pvals) in enumerate(pval_items[1:]) :
            b = bar_ax.bar(bins[:-1],pvals,bottom=sum_bottoms,width=bar_width,color=colors[i+1])
            bars.append(b)
            sum_bottoms += pvals
        bar_ax.legend([b[0] for b in bars],[x[0] for x in pval_items])
        bar_ax.axis((-10,max(bins),0,max(sum_bottoms)))
        bar_ax.set_title('Peak map distribution by pvalue')
        bar_ax.set_xlabel('-10*log10(pvalue)')
        bar_ax.set_ylabel('relative log10(# peaks)')

        pval_fn = fn_base+'_pval_bar.png' if opts.bar_fn is None else opts.bar_fn
        mp.savefig(pval_fn)
