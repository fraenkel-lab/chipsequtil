
import random
import re
import sys
from collections import defaultdict

from chipsequtil import get_org_settings, get_gc_content, get_gc_content_distribution, RefGeneFile
from nib import NibDB, NibException
from TAMO.seq import Fasta

def rejection_sample_bg(fg_dict,organism,bins=100,num_samples=None,verbose=False) :
    '''Generate background sequences according to the size, distance from genes,
    and GC content distributions of the supplied foreground sequences.  *fg_dict*
    is a dictionary of <header>:<sequence> items, where the first part of the
    header must contain:

    >chrX:<start>-<end>

    *organism* is a string that will be used to call the *chipsequtil.get_org
    settings* function and uses the 'genome_dir' and 'annotation_path' keys.
    *bins* is the number of bins to use for representing the GC content
    distribution.  Function returns a dictionary of <header>:<sequence> items
    of generated background sequences.'''

    nib_db = NibDB(nib_dirs=[get_org_settings(organism)['genome_dir']])
    tss_fn = get_org_settings(organism)['annotation_path']
    tss = defaultdict(list)
    for rec in RefGeneFile(tss_fn) :
        tss[rec['chrom']].append((int(rec['txStart']),int(rec['txEnd']),))

    # for each peak find the chromosome, distance to nearest
    # gene, size of peaks in bases, and GC content
    num_samples = len(fg_dict) if not num_samples else num_samples
    dists,sizes=[],[]

    for header,seq in fg_dict.items() :

        # chromosome first field in fasta headers from bed2seq.bedtoseq
        chrom = header.split(':')[0]

        # adjust chromosomes in special cases
        if re.search('random',chrom.lower()) or chrom.lower() == 'chrm' :
            continue
        if chrom == 'chr20' :
            chrom = 'chrX'
        elif chrom == 'chr21' :
            chrom = 'chrY'

        # start first int in second field of bed2seq.bedtoseq header
        start = int(header.split(':')[1].split('-')[0])
        midpoint = start + len(seq)/2

        # figure out which chromosome we're working on
        tss_chr = tss[chrom]

        # dsts_to_genes is the distance of this peak from all the genes, find minimum
        dists_to_genes = [(s[0]-midpoint) for s in tss_chr]
        try :
            min_dist = min(dists_to_genes,key=lambda x : abs(x))
        except :
            err_str = '''Warning: no genes were found for sequence with header
                         %s, not using to calculate distributions.
                         '''%header
            sys.stderr.write(err_str)
        dists.append(min_dist)

        # calculate # bases
        sizes.append(len(seq))

    # GC content distribution for the foreground sequences
    gc_dist = get_gc_content_distribution(fg_dict.values(),bins=bins)

    # max_gc is # peaks w/ highest GC content
    max_gc = max(gc_dist)

    # start generating bg sequences
    bg_dict = {}

    bg_gcs,bg_sizes=[],[]

    # gene_starts is a list of all genes in (chromosome,gene start) tuples
    gene_starts=[]
    for key in tss.keys():
        chrom=key.split('chr')[-1]
        for x in tss[key]:
            gene_starts.append((key,x[0]))

    # generate a bg sequence for every fg sequence
    for i in range(num_samples):
        if verbose : sys.stderr.write('\n%d/%d'%(i,num_samples))

        # propose sequences until one is accepted
        accepted_sequence = False
        while not accepted_sequence:
            if verbose : sys.stderr.write('. ')

            # sample a random distance from the list of distances
            d = random.choice(dists)

            # pick a random gene
            chrom, coord = random.choice(gene_starts)

            # propose a starting point for the bg sequence
            midpoint = coord-d+random.randint(-100,100)

            # propose a size for the bg sequence
            size = random.choice(sizes)
            start = int(midpoint-int(size/2))
            stop = int(midpoint+int(size/2))

            # if start or stop are negative, skip and try again
            if start < 0 or stop < 0 : continue

            # randomly choose strand
            strand = '+' if random.random() > 0.5 else '-'

            # extract the proposed sequence
            try :
                nib_title, seq = nib_db.get_fasta(chrom,start,stop,strand)
            except IOError :
                if verbose : sys.stderr.write('IOError in NibDB, skipping: %s,%d-%d,%s\n'%(chrom,start,stop,strand))
                continue
            except NibException :
                if verbose : sys.stderr.write('NibDB.get_fasta error\n')
                continue

            # determine the GC bin for this sequence
            gc_content = get_gc_content(seq)
            gc_bin = -1
            for i in range(bins) :
                win_start = i/float(bins)
                win_end = (i+1)/float(bins)
                if gc_content >= win_start and gc_content < win_end :
                    gc_bin = i
                    continue

            # pick a uniform random number such that it does not exceed
            # the maximum GC content distribution over bins
            # if the random number is <= the GC content for this
            # proposed sequence, accept, otherwise reject
            r = random.random() * max_gc
            if r > gc_dist[gc_bin] :
                continue
            else:
                bg_gcs.append(x)
                bg_sizes.append(size)
                accepted_sequence = True
                header = '%s:%d-%d'%(chrom,start,stop)
                bg_dict[header] = seq

    return bg_dict
