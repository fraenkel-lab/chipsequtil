#!/usr/bin/env python

import math
import sys
from csv import reader, DictWriter
from optparse import OptionParser
from scipy.stats import zprob, std

from TAMO2.MotifTools import load

usage = "%prog [options] <fg motif filename> <THEME randomization filename>"
description = """
Calculate adjusted p-values for enrichmed motifs found by THEME2.py.  <fg motif filename>
is the output file corresponding to the --motif-filename option of THEME2. <THEME
randomization file> is the output file corresponding to the --random-output
THEME2 option (when THEME2 was invoked with --randomization flag). p-values are
adjusted using Benjamini-Hochberg.
"""

parser = OptionParser(usage=usage,description=description)
parser.add_option('--output',dest='output',default=None,help='filename to write output to')


def benjamini_hochberg_adjust(pvals) :
    """Perform Benjamini Hochberg multiple testing pvalue adjustment.  Returns
    list of adjusted pvalues in the same order as they were passed."""

    # need to sort pvalues descending and maintain original ordering
    # to reconstruct later
    pvals_ranks = zip(pvals,range(len(pvals)))
    pvals_ranks.sort(reverse=True)
    orig_pvals, orig_ranks = zip(*pvals_ranks)

    n_pvals = float(len(pvals))

    adj_pvals = []
    cum_min = 1.0
    for i, pval in enumerate(orig_pvals) :
        rank_factor = n_pvals/(n_pvals-i)

        # we take the min of the smallest seen and adjusted pvalue
        # because we interpret BH p-values as FDR values, thus
        # the we have confidence in the lowest FDR we find regardless
        # of rank
        adj_pval = min(cum_min, pval*rank_factor)
        adj_pvals.append(adj_pval)

    # put the pvalues in their original order
    adj_pvals_ranks = zip(orig_ranks,adj_pvals)
    adj_pvals_ranks.sort()

    ranks, adj_pvals = zip(*adj_pvals_ranks)

    return adj_pvals


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 2 :
        parser.error('Must provide exactly two non-option arguments')

    fg_motifs_fn, random_motifs_fn = args

    fg_motifs = load(fg_motifs_fn)

    random_motifs = []
    for rec in reader(open(random_motifs_fn),delimiter='\t') :
        random_motifs.append((rec[0],float(rec[1]),eval(rec[2])))

    # calculate statistics on motifs
    assert len(fg_motifs) == len(random_motifs)
    zipped_motifs = zip(fg_motifs,random_motifs)

    motifs = []
    pvals = [] # for later multiple testing adjument
    motif_headers = ['motif_name','motif','pvalue','cverror','mean_rand_cverror','rand_stddev','beta','threshold']
    for motif, rand_stats in zipped_motifs :
        motif_str, mean_cv, samples = rand_stats

        # calculate standard deviation
        rand_stddev = std(samples)

        # calculate zscore and invert to pvalue
        zscore = (motif.cverror - mean_cv)/rand_stddev
        pval = zprob(zscore)

        try :
            # motif 'name' is buried in Motif.source
            id, motif_name, some_other_motif_name, rest = motif.source.split('\t',3)
        except ValueError :
            # probably the source attribute doesn't have the last field
            id, motif_name, some_other_motif_name = motif.source.split('\t',3)

        motif_dict = dict(zip(motif_headers,[some_other_motif_name,motif,pval,motif.cverror,mean_cv,rand_stddev,motif.beta,motif.match_thresh]))
        motifs.append(motif_dict)
        pvals.append(pval)

    # do multiple testing adjustment
    adj_pvals = benjamini_hochberg_adjust(pvals)
    for s, adj_pval in zip(motifs, adj_pvals) :
        s['pvalue'] = adj_pval

    motifs.sort(key=lambda x: x['pvalue'])

    # write out results
    out_f = open(opts.output,'w') if opts.output else sys.stdout
    motif_writer = DictWriter(out_f,delimiter='\t',fieldnames=motif_headers,lineterminator='\n')
    motif_writer.writerow(dict(zip(motif_headers,motif_headers)))
    motif_writer.writerows(motifs)
    if opts.output : out_f.close()
