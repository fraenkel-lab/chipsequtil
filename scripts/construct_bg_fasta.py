#!/usr/bin/env python

import os
import sys
import warnings

from collections import defaultdict
from optparse import OptionParser

from chipsequtil import get_org_settings, RefGeneFile
from chipsequtil.nib import NibDB
from chipsequtil.util import MultiLineHelpFormatter
from TAMO.seq import Fasta

usage='%prog [options] <type> <organism> <foreground fasta>'
description='Create background sequence databses for motif finding, etc.'
parser = OptionParser(usage=usage,description=description,formatter=MultiLineHelpFormatter())


def rejection_sampling(fg,settings_dict,gc_bins=20) :

    genm_db = NibDB(settings_dict['genome_dir'])
    annot = RefGeneFile(settings_dict['annotation_file'])


    num_peak_bases = 0
    for header, seq in fg.items() :
        num_peak_bases += len(seq)


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 3 :
        parser.error('Must provide three non-option arguments')

    sample_type, organism, fg_fn = args[:3]

    settings_dict = get_org_settings(organism)

    fg = Fasta.load(fg_fn)
    bg = rejection_sampling(fg,settings_dict)


###############################################################
# start Chris' code from rej_samp_bg_rand2.py
    the_genes={} #list of distances to nearest TSS

    # for each peak find the chromosome, distance to nearest
    # gene, size of peaks in bases, and GC content
    the_chrs,dists,sizes,gcs=[],[],[],[]

    # number of bases in the fg sequences
    size=0

    for key in pos_seqs.keys():

        size+=len(pos_seqs[key])

        # chromosome first field in fasta headers from bed2seq.bedtoseq
        chr=key.split(':')[0]

        # adjust chromosomes in special cases
        if re.search('random',chr):
            continue
        if chr=='chr20':
            chr='chrX'
        elif chr=='chr21':
            chr='chrY'
        if not the_genes.has_key(chr):
            the_genes[chr]=[]

        # start first int in second field of bed2seq.bedtoseq header
        start=int(key.split(':')[1].split('-')[0])
        midpoint=int(start+len(pos_seqs[key])/2)

        # figure out which chromosome we're working on
        tss_chr=tss[chr.split('chr')[-1]]

        # D is the distances from all the genes, find minimum
        D=[(s[0]-midpoint) for s in tss_chr]

        # best distance for this peak
        minD=min([abs(x) for x in D])
        best=[d for d in D if abs(d)==minD]
        dists.append(best[0])

        # chromosome for this peak
        the_chrs.append(chr)
        seq=pos_seqs[key]

        # calculate # bases and GC content
        N=len(seq)
        sizes.append(N)
        gc=len([x for x in seq if (x=='G')or(x=='C')])/N
        gcs.append(gc)

    #bin GC content distribution
    bins=20

    # q is # of peaks w/ x% GC content
    q=[0]*bins

    for gc in gcs:
        for i in range(bins):
            win_start=i/bins
            win_end=(i+1)/bins
            if gc>=win_start and gc<win_end:
                q[i]+=1
                continue

    # q is now % peaks w/ x% GC content
    q=[x/Nseqs for x in q]
    #print q

    # c is # peaks w/ highest GC content
    c=max(q)*Nseqs

    # start generating bg sequences
    print "Done assembling distance and gc content distributions"
    genome_outfile=open(bg,'w')

    # make twice as many 
    size=round(size/(2*len(pos_seqs)))
    bg_gcs,bg_sizes=[],[]
    #for key in the_genes.keys():
        #chrom=key.split('chr')[-1]
        #the_genes[key]=[x[0] for x in tss[chrom]]

    # C_TX is a list of all genes in (chromosome,gene start) tuples
    C_TX=[]
    for key in tss.keys():
        chrom=key.split('chr')[-1]
        for x in tss[chrom]:
            C_TX.append((chrom,x[0]))

    # generate a bg sequence for every fg sequence
    for i in range(Nseqs):

        # propose sequences until one is accepted
        keep_going=1
        while keep_going:
            #random.shuffle(the_chrs)

            # randomize the list of distances from genes
            random.shuffle(dists)
            #chr=the_chrs[0]

            # pick the first distance, i.e. at random
            d=dists[0]

            #random.shuffle(the_genes[chr])

            # randomize the gene list
            random.shuffle(C_TX)

            # randomize the peak sizes
            random.shuffle(sizes)

            # pick a random gene
            (chr,coord)=C_TX[0]

            #coord=the_genes[chr][0]
            # propose a starting point for the bg sequence
            midpoint=coord-d+random.randint(-100,100)

            # propose a starting size for the bg sequence
            size=sizes[0]
            start=int(midpoint-int(size/2))
            stop=int(midpoint+int(size/2))
            id='chr'+chr.split('chr')[-1]+':'+str(start)+'-'+str(stop)
            r=random.random()

            # randomly choose strand
            if r<0.5: strand='+'
            else: strand='-'

            # extract the proposed sequence
            nib_title,seq=nibfrag.sequence('chr'+chr,start, stop,strand)
            if not seq:
                print 'NOT FOUND', chr,start,stop,
                continue
            else:

                N,y=0,0
                # calculate the GC content for the proposed sequence
                for line in seq:
                    s=line.upper()
                    N+=len(line)
                    y+=len([x for x in s if (x=='G')or(x=='C')])
                    if line[0]=='N': continue
                x=float(y)/N

                # determine the GC bin for this sequence
                #gc=float(len([x for x in seq if (x=='G')or(x=='C')]))/N
                for i in range(bins):
                    win_start=i/bins
                    win_end=(i+1)/bins
                    if x>=win_start and x<win_end:
                        bin=i
                        continue

                # pick a uniform random number such that it does not exceed
                # the maximum GC content distribution over bins
                r=random.random()*c/Nseqs

                # if the random number is <= the GC content for this
                # proposed sequence, accept, otherwise reject
                if r>q[bin]:
                    #print 'skip'
                    continue
                else:
                    #print bin
                    bg_gcs.append(x)
                    bg_sizes.append(size)
                    keep_going-=1
                    title='>%s\n'%id
                    genome_outfile.write(title)
                    for line in seq:
                        genome_outfile.write(line.upper()+'\n')
    print len(gcs)
    print len(bg_gcs)
    fg_mean,fg_sdev=mean_sdev(gcs)
    print fg_mean,fg_sdev
    #bg_mean,bg_sdev=mean_sdev(bg_gcs)
    bg_mean=scipy.mean(bg_gcs)
    bg_sdev=scipy.std(bg_gcs)
    print bg_mean,bg_sdev
    fg_size_m,fg_size_dev=mean_sdev(sizes)
    bg_size_m,bg_size_dev=mean_sdev(bg_sizes)
    print fg_size_m,fg_size_dev
    print bg_size_m,bg_size_dev
    genome_outfile.close()

