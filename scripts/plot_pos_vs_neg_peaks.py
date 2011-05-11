#!/usr/bin/env python

import os
import sys

import matplotlib
matplotlib.use('AGG')

from matplotlib.pyplot import *
from numpy import arange, log10
from optparse import OptionParser

from chipsequtil import MACSFile

usage = '%prog [options] <pos peaks fn> <neg peaks fn>'
parser = OptionParser(usage=usage)
parser.add_option('-o','--output',dest='out_fn',default=None,help='filename of output image')

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])
    pos_fn, neg_fn = args

    pos_f, neg_f = MACSFile(pos_fn), MACSFile(neg_fn)

    pos_peaks = []
    pos_pvals = []
    for pk in pos_f :
        pos_pvals.append(float(pk['-10*log10(pvalue)'])/10.)
        pos_peaks.append((pk['-10*log10(pvalue)'],pk))

    pos_peaks.sort()

    neg_peaks = []
    neg_pvals = []
    for pk in neg_f :
        neg_pvals.append(float(pk['-10*log10(pvalue)'])/10.)
        neg_peaks.append((pk['-10*log10(pvalue)'],pk))

    neg_peaks.sort()

    min_pval, max_pval = min(pos_pvals+neg_pvals), max(pos_pvals+neg_pvals)

    pval_rng = arange(min_pval,max_pval,(max_pval-min_pval)/100.)

    # construct cdfs
    pos_cdf, neg_cdf = [], []
    for pval in pval_rng :
        pos_cdf.append(len(filter(lambda x: x >= pval,pos_pvals)))
        neg_cdf.append(len(filter(lambda x: x >= pval,neg_pvals)))

    # normalize cdfs
    pos_cdf_norm = [1.*x/max(pos_cdf) for x in pos_cdf]
    neg_cdf_norm = [1.*x/max(neg_cdf) for x in neg_cdf]

    # log of pvals
    pos_logs = map(log10,pos_cdf)
    neg_logs = map(log10,neg_cdf)
    plot(pval_rng,pos_logs)
    plot(pval_rng,neg_logs)
    ytics, ylabs = yticks()
    clf()

    # normalize logs for plotting
    pos_logs_norm = [1-1.*x/max(pos_logs) for x in pos_logs]
    neg_logs_norm = [1-1.*x/max(neg_logs) for x in neg_logs]

    # calculate pos proportion for each pvalue
    pos_ratio = []
    pos_only = []
    for pos, neg in zip(pos_cdf_norm,neg_cdf_norm) :
        #pos_ratio.append(pos/(pos+neg))
        if neg == 0 :
            pos_only.append(pos_ratio[-1])
            #pos_ratio.append(pos_ratio[-1])
        else :
            pos_ratio.append(pos/neg)

    subplot(211)
    plot(pval_rng, pos_logs, 'b-')
    plot(pval_rng, neg_logs, 'g-')
    yticks(ytics,[int(10**y) for y in ytics])
    title('positive vs. negative peaks')
    legend(('positive','negative'),loc='upper right')
    xlabel('-log(p-value)')
    ylabel('# Peaks')
    axis('tight')

    subplot(212)
    plot(pval_rng[:len(pos_ratio)], map(log10,pos_ratio), 'k-')
    plot(pval_rng[len(pos_ratio):], map(log10,pos_only),'k--')
    #plot(pval_rng,pos_ratio, 'k-')
    axis('tight')
    xlabel('-log(p-value)')
    #ylabel('# pos / (# pos + # neg)')
    ylabel('log(# pos / # neg)')

    if opts.out_fn is None :
        pos_base_fn, pos_fn_ext = os.path.splitext(pos_fn)
        out_fn = '%s_pos_v_neg.png'%pos_base_fn
    else :
        out_fn = opts.out_fn
    savefig(out_fn)
