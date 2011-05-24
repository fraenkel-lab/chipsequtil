#!/usr/bin/env python

import matplotlib
matplotlib.use('AGG')

import numpy as np
import os
import random
import string
import sys

from math import log, pow
import matplotlib.pyplot as mp
from multiprocessing import Pool
from optparse import OptionParser 
from scipy.stats.stats import pearsonr

from chipsequtil import MACSFile, get_org_settings
from chipsequtil.nib import NibDB
from chipsequtil.sampling import rejection_sample_bg
from TAMO import MotifTools as mt
from TAMO.MotifTools import load

usage = "%prog [options] <org> <peaks fn> <TAMO motif fn>"
desc = "Do some motif scanning stuffs"
parser = OptionParser(usage=usage,description=desc)

parser.add_option('-n','--top-n',dest='top_n',type='int',default=None,help='use top n peaks by pvalue for sequence scanning [default: all]')
parser.add_option('-i','--motif-indices',dest='motif_ind',default='all',help='which indices from <TAMO motif fn> to use [default: %default]')
parser.add_option('-d','--dir',dest='dir',default='motif_results',help='write all results into this directory')
parser.add_option('--fixed-peak-width',dest='fixed_w',type='int',default=None,help='use only a fixed peak window around the summit instead of whole peak')

revcomp_map = string.maketrans('ACGT','TGCA')

def score_sequence(seq,motif) :
    ll_max = -sys.maxint
    for i in range(len(seq)-len(motif)) :
        # forward strand
        ll_for_sum = 0
        subseq = seq[i:i+len(motif)].upper()
        for n,pos in zip(subseq,motif.ll) :
            ll_for_sum += pos[n]
        # reverse strand
        ll_rev_sum = 0
        subseq = reversed(subseq.translate(revcomp_map))
        for n,pos in zip(subseq,motif.ll) :
            ll_rev_sum += pos[n]
        ll_max = max(ll_max,ll_for_sum,ll_rev_sum)

    return ll_max

illegal_fn_chars = '/;& ()'
fn_trans = string.maketrans(illegal_fn_chars,'_'*len(illegal_fn_chars))

def fasta_itr(fn) :
    f = open(fn)
    header = None
    seq = None
    for l in f :
        if l.strip().startswith('>') :
            if seq is not None :
                yield (header,seq)
                seq = None
            header = l.strip()
        else :
            seq = seq+l.strip() if seq is not None else l.strip()

    # last record
    yield (header, seq)

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 3 :
        parser.error('Exactly 3 non-option arguments must be provided')

    org, peaks_fn, motif_fn = args

    if not os.path.exists(opts.dir) :
        os.mkdir(opts.dir)

    peaks_dt = np.dtype([('chr',np.str_,13),('start',np.int32),('end',np.int32),('pvalue',np.float64)])
    if opts.fixed_w is not None :
        
        peaks = np.array([(r['chr'],
                          r['start']+r['summit']-opts.fixed_w/2.,
                          r['start']+r['summit']+opts.fixed_w/2.,
                          r['-10*log10(pvalue)']) for r in MACSFile(peaks_fn)],
                          dtype=peaks_dt)
    else :
        peaks = np.array([(r['chr'],
                           r['start'],
                           r['end'],
                           r['-10*log10(pvalue)']) for r in MACSFile(peaks_fn)],
                           dtype=peaks_dt)

    # -10*log10(pvalue) -> -log10(pvalue)
    peaks[:]['pvalue'] /= 10.
    peak_pvals = peaks[:]['pvalue']

    # find the sorted order of peaks by descending pvalue
    peak_pval_inds = peak_pvals.argsort()
    peak_pval_inds = peak_pval_inds[::-1] # ascending -> descending
    peaks = peaks[peak_pval_inds,:]

    if opts.top_n is not None :
        peaks = peaks[0:opts.top_n]
        peak_pvals = peak_pvals[peak_pval_inds][0:opts.top_n]

    # extract fasta sequences for these peaks
    nibDb = NibDB(nib_dirs=get_org_settings(org)['genome_dir'])

    # get the peak sequences
    sys.stderr.write('Getting peak sequences\n')
    fasta_batch = []
    for i in range(peaks.size) :
        fasta_batch.append((str(peaks[i]['chr']),int(peaks[i]['start']),int(peaks[i]['end']),'+'))
    fg_fasta_headers, fg_fasta = nibDb.get_fasta_batch(fasta_batch)

    # need a dict for background sampling
    # headers have genome_dir and .nib in them, strip that out
    sys.stderr.write('Converting nib output to dict\n')
    fg_fasta_headers = list(fg_fasta_headers)
    fg_fasta_dict = {}
    for h,s in zip(fg_fasta_headers,fg_fasta) :
        h = h.replace('>'+get_org_settings(org)['genome_dir']+'/','')
        h = h.replace('.nib','')
        if len(s) > 150 :
            fg_fasta_dict[h] = s

    # now sample the background sequences
    sys.stderr.write('Sampling bg sequences (len(fg_fasta)==%d)\n'%(len(fg_fasta_dict)))
    bg_fasta_dict = rejection_sample_bg(fg_fasta_dict,org,bg_match_epsilon=1e-3,verbose=True)
    bg_fasta = bg_fasta_dict.values()

    # load the motifs
    sys.stderr.write('Movin right along\n')
    motifs = load(motif_fn)

    if opts.motif_ind != 'all' :
        motif_indices = [int(i) for i in opts.motif_ind.split(',') if len(i) != 0]
        motifs = [motifs[i] for i in motif_indices]
    else :
        motif_indices = range(len(motifs))

    # use all cores w/ a Pool
    #pool = Pool(processes=opts.n_procs)

    # go through each motif
    job_params = []
    res = []
    #for i,m in zip(motif_indices,motifs) :
    #    job_params.append((i,m,peak_pvals,fg_fasta,bg_fasta,opts.dir))
    #seq_scores = pool.map(analyze_motif_sequences,job_params)

    seq_scores = []
    for m_i,m in zip(motif_indices,motifs) :
        out_dir = opts.dir

        try :
            m_name = m.source.split('\t')[2]
        except :
            m_name = m.source.split()[0]

        print 'starting',m_name

        fg_ratios = []
        for seq in fg_fasta :
            #max_score = score_sequence(seq,m)
            max_score = m.bestscan(seq.upper())
            fg_ratios.append((max_score-m.minscore)/(m.maxscore-m.minscore))
        fg_ratios = np.array(fg_ratios)

        bg_ratios = []
        for seq in bg_fasta :
            #max_score = score_sequence(seq,m)
            max_score = m.bestscan(seq.upper())
            bg_ratios.append((max_score-m.minscore)/(m.maxscore-m.minscore))
        bg_ratios = np.array(bg_ratios)

        fg_mean = sum(fg_ratios)/len(fg_ratios)
        fg_std = np.std(fg_ratios)
        bg_mean = sum(bg_ratios)/len(bg_ratios)
        bg_std = np.std(bg_ratios)

        m_mat = np.array((fg_ratios,bg_ratios,peak_pvals))
        fg_score_sort_inds = m_mat[0,:].argsort()

        motif_score_cnts, motif_score_bins = np.histogram(m_mat[0,:],bins=20)
        binned_motif_scores = []
        for st, end in zip(motif_score_bins[:-1],motif_score_bins[1:]) :
            binned_motif_scores.append(m_mat[2,(m_mat[0,:]>=st)&(m_mat[0,:]<end)])

        mp.figure(figsize=(4,4))
        font = {'size':'9'}
        mp.rc('font',**font)

        mp.plot(fg_ratios,peak_pvals,'bo')

        # calculate pearson correlation coefficient
        pear_r, pear_pval = pearsonr(fg_ratios,peak_pvals)
        mp.title('Max motif strength vs peak pvalue\n(r=%.2f,pval=%.2g)'%(pear_r,pear_pval))
        img_fn = os.path.join(out_dir,m_name.translate(fn_trans)+'_%d_corr.png'%m_i)
        mp.savefig(img_fn)
        mp.clf()

        # line plot of average peak p-value for binned motif score
        mp.title('Average peak p-value for binned motif score\n%s'%m_name)
        mp.xlabel('normalized motif score')
        mp.ylabel('avg(-log10(pvalue))')
        mp.boxplot(binned_motif_scores,positions=np.arange(motif_score_bins.size-1),sym='')
        p = mp.plot(np.arange(motif_score_bins.size-1),
                [x.mean() for x in binned_motif_scores],
                'bo',
                label='Mean fg score')
        p = p[0]

        # draw a crosshair
        bg_median_ind = np.argwhere(((motif_score_bins<=bg_mean)[:-1] & (motif_score_bins>=bg_mean)[1:])).ravel()[0]
        bg_median = np.median(binned_motif_scores[bg_median_ind])
        xlim, ylim = p.axes.get_xlim(), p.axes.get_ylim()
        mp.plot([bg_median_ind,bg_median_ind],ylim,'k-',label='Mean bg score=%.2g'%m_mat[1,:].mean())
        mp.plot(xlim,[bg_median,bg_median],'k-')
        mp.xticks(np.arange(motif_score_bins.size)[1::5],['%.2f'%x for x in motif_score_bins[1::5]])
        mp.legend(loc='upper left')

        img_fn = os.path.join(out_dir,m_name.translate(fn_trans)+'_%d_peakmot.png'%m_i)
        mp.savefig(img_fn)
        mp.clf()

        ret_d ={'m_name': m_name,
                'fg_mean': fg_mean,
                'fg_std': fg_std,
                'bg_mean': bg_mean,
                'bg_std': bg_std,
                'fg_scores': fg_ratios,
                'bg_scores': bg_ratios,
                #'wmw_pval': WMWtest(fg_ratios,bg_ratios)
               }

        print 'done with',m_name

        seq_scores.append(ret_d)
