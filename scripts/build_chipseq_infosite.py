#!/usr/bin/env python

import getpass
import glob
import json
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as mp
import os
import re
import shutil
import sys

from collections import defaultdict
from csv import reader, writer, DictReader
from math import log
from optparse import OptionParser
from subprocess import call

from chipsequtil import MACSFile, get_org_settings
from reStUtil import *

usage = '%prog [options] [<peak filename> <peak filename> ...]'
parser = OptionParser(usage=usage)
parser.add_option('-d','--dir',dest='dir',default='.',help='Source directory [default: %default]')
parser.add_option('-n','--name',dest='name',help='Experiment name [default: current directory name]')
parser.add_option('--skip-motif-scan',dest='skip_motif_scan',action='store_true',help="skip motif_scan.py, but still build motifs into document (assumes motif_scan.py was previously run)")
parser.add_option('--skip-motif-stuff',dest='skip_motif_stuff',action='store_true',help="motif stuff takes a long time, manually skip it if no motif results are available or you don't care about them")

{
  "experiment path": "/nfs/antdata/analysis/100809_P/100809_St7-10ul-p300_degenhar_Fraenkel_L2_mm9_sorted.bed", 
  "analysis path": "/net/ventral/nfs/people/labadorf/analysis/100809_P_St7_10ul", 
  "stage url": "http://fraenkel.mit.edu/stage/labadorf", 
  "peak files": {
    "100809_P_St7_10ul_mfold10,30_pval1e-5": {
      "total tags in control": 9331149, 
      "total tags in treatment": 10064908, 
      "Range for calculating regional lambda": "1000 bps and 10000 bps", 
      "tag size": 35, 
      "name": "100809_P_St7_10ul_mfold10,30_pval1e-5", 
      "model fold": "10,30", 
      "format": "BED", 
      "tags after filtering in treatment": 5099883, 
      "band width": 150, 
      "Redundant rate in control": 0.40999999999999998, 
      "Redundant rate in treatment": 0.48999999999999999, 
      "effective genome size": 2110000000.0, 
      "d": 145, 
      "maximum duplicate tags at the same position in control": 1, 
      "control file": "cntrl_6-3_sorted_filterbed.txt", 
      "MACS version": "1.4.0beta", 
      "ChIP-seq file": "exp_100809_St7-10ul-p300_degenhar_Fraenkel_L2_mm9_sorted.bed", 
      "tags after filtering in control": 5481613, 
      "maximum duplicate tags at the same position in treatment": 2, 
      "pvalue cutoff": 1.0000000000000001e-05
    }
  }, 
  "format": "BED", 
  "FDR filter": "none", 
  "experiment name": "100809_P_St7_10ul", 
  "mapping type": "TSS", 
  "pipeline args": {
    "--filter-peaks-args": "--sort-by=pvalue --top=200", 
    "--macs-args": "--mfold=10,30 --tsize=35 --bw=150 --format=BED --pvalue=1e-5", 
    "--map-args": "--tss --upstream-window=10000 --downstream-window=10000"
  }, 
  "org": "mm9", 
  "control path": "/nfs/antdata/analysis/090828_42JVC/6-3/6-3_sorted_filterbed.txt", 
  "mapping window": [
    "10000", 
    "10000"
  ], 
  "peaks used by THEME": "200", 
  "stage_dir": "/nfs/antdata/web_stage/labadorf"
}

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    exp_dir = os.path.abspath(opts.dir)
    exp_name = opts.name if opts.name is not None else os.path.basename(exp_dir)

    # 1. find the param JSON file
    param_json_fn = glob.glob('*params.json')
    if len(param_json_fn) == 0 :
        sys.stderr.write('Could not find parameter file, building one as best I can\n')
        curr_user = getpass.getuser()
        json_d = {'analysis path':os.getcwd(),
                  'stage url':'http://fraenkel.mit.edu/stage/'+curr_user,
                  'stage dir':'/nfs/antdata/web_stage/'+curr_user
                 }
    else :
        if len(param_json_fn) > 1 :
            sys.stderr.write('Found more than one parameter file, picking the first one: %s\n'%','.join(param_json_fn))
        param_json_fn = param_json_fn[0]
        json_d = json.load(open(param_json_fn))

    # 2. make a new directory to save all the stuff
    infosite_dir_name = exp_name+'_infosite'
    infosite_path = os.path.join(os.getcwd(),infosite_dir_name)
    if not os.path.exists(infosite_path) :
        os.mkdir(infosite_path)

    infosite_img_path = os.path.join(infosite_path,'images')
    if not os.path.exists(infosite_img_path) :
        os.mkdir(infosite_img_path)

    # 3. setup web staging directory
    stage_dir_path = os.path.join(json_d['stage dir'],infosite_dir_name)
    if not os.path.exists(stage_dir_path) :
        os.symlink(infosite_path,stage_dir_path)

    # 4. get the peaks files stats, don't want negative peaks
    if len(args) == 0 :
        peaks_fns = glob.glob('*_peaks.xls')
        peaks_fns = filter(lambda x: 'negative' not in x,peaks_fns)
    else :
        peaks_fns = args
    analysis_sets = []
    peak_json = json_d['peak files'] = {}

    # analyze all the peak files
    for peak_fn in peaks_fns :
        print 'processing:',peak_fn
        macs_f = MACSFile(peak_fn)
        peak_json[peak_fn] = macs_f.file_info

        # positive peaks
        peak_stats = defaultdict(list)
        num_peaks = 0
        pos_chr_dist = defaultdict(int)
        for peak in macs_f :
            pos_chr_dist[peak['chr']] += 1
            peak_stats['length'].append(peak['length'])
            peak_stats['tags'].append(peak['tags'])
            peak_stats['pvalue'].append(peak['-10*log10(pvalue)'])
            peak_stats['fold_enrichment'].append(peak['fold_enrichment'])
            peak_stats['fdr'].append(peak['FDR(%)'])
            num_peaks += 1

        peak_json[peak_fn]['positive peaks'] = num_peaks
        peak_json[peak_fn]['reads under peaks'] = sum(peak_stats['tags'])

        # extract paired peaks info out of output.txt
        output_fn = peak_json[peak_fn]['name']+'_output.txt'
        output_regexes = ('#2 number of (paired peaks): (\d+)',)
        for l in open(output_fn) :
            for regex in output_regexes :
                m = re.search(regex,l)
                if m is not None :
                    peak_json[peak_fn][m.group(1)] = int(m.group(2))

        # do the negative peaks
        # negative peak file is now filtered
        neg_peak_fns = glob.glob(peak_json[peak_fn]['name']+'_negative_peaks_*.xls')
        neg_peak_stats = None

        #TODO - do check for file exists
        if neg_peak_fns :
            neg_peak_fn = neg_peak_fns[0]
            neg_peak_f = MACSFile(neg_peak_fn)

            neg_peak_stats = defaultdict(list)
            num_peaks = 0
            neg_chr_dist = defaultdict(int)
            for peak in neg_peak_f :
                neg_chr_dist[peak['chr']] += 1
                neg_peak_stats['length'].append(peak['length'])
                neg_peak_stats['tags'].append(peak['tags'])
                neg_peak_stats['pvalue'].append(peak['-10*log10(pvalue)'])
                neg_peak_stats['fold_enrichment'].append(peak['fold_enrichment'])
                neg_peak_stats['fdr'].append(peak['FDR(%)'])
                num_peaks += 1

            peak_json[peak_fn]['negative peaks'] = num_peaks
            peak_json[peak_fn]['reads under negative peaks'] = sum(peak_stats['tags'])
        else :
            peak_json[peak_fn]['negative peaks'] = 'NA'
            peak_json[peak_fn]['reads under negative peaks'] = 'NA' 

        # save the track lines
        ucsc_track_fn = peak_json[peak_fn]['name']+'_MACS_wiggle_tracks.txt'
        if os.path.exists(ucsc_track_fn) :
            peak_json[peak_fn]['ucsc tracks'] = open(ucsc_track_fn).readlines()

        #font = {'size':'9'}
        #mp.rc('font',**font)

        figsize = (3.5,3.5)
        subplots_sizes = {'top':0.8,'left':0.15,'right':0.95}

        # sometimes there are no negative peaks
        hist_labels = ['+ peaks']
        hist_peak_stats = [peak_stats]
        if neg_peak_stats is not None :
            hist_labels.append('- peaks')
            hist_peak_stats.append(neg_peak_stats)

        #hist_labels = ('+ peaks','- peaks')
        # create histograms for each of the attributes
        len_hist_name = macs_f.file_info['name']+'_length.png'
        len_hist_fn = os.path.join(infosite_img_path,len_hist_name)
        len_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+len_hist_name
        peak_json[peak_fn]['length distribution url'] = len_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)

        hist_stats = [p['length'] for p in hist_peak_stats]
        mp.hist(hist_stats,label=hist_labels,bins=20,log=True)
        mp.title('%s\npeak length distribution'%macs_f.file_info['name'])
        mp.xlabel('peak length')
        mp.ylabel('# peaks')
        mp.legend()
        print len_hist_fn
        mp.savefig(len_hist_fn)
        mp.clf()

        tags_hist_name = macs_f.file_info['name']+'_tags.png'
        tags_hist_fn = os.path.join(infosite_img_path,tags_hist_name)
        tags_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+tags_hist_name
        peak_json[peak_fn]['tag distribution url'] = tags_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)
        hist_stats = [p['tags'] for p in hist_peak_stats]
        mp.hist(hist_stats,label=hist_labels,bins=20,log=True)
        mp.title('%s\npeak tag count distribution'%macs_f.file_info['name'])
        mp.xlabel('# tags')
        mp.ylabel('# peaks')
        mp.legend()
        mp.savefig(tags_hist_fn)
        mp.clf()

        pval_hist_name = macs_f.file_info['name']+'_pval.png'
        pval_hist_fn = os.path.join(infosite_img_path,pval_hist_name)
        pval_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+pval_hist_name
        peak_json[peak_fn]['pvalue distribution url'] = pval_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)
        hist_stats = [p['pvalue'] for p in hist_peak_stats]
        mp.hist(hist_stats,label=hist_labels,bins=20,log=True)
        mp.title('%s\npeak -10*log10(p-valuek) distribution'%macs_f.file_info['name'])
        mp.xlabel('-10*log10(p-value)')
        mp.ylabel('# peaks')
        mp.legend()
        mp.savefig(pval_hist_fn)
        mp.clf()

        fold_hist_name = macs_f.file_info['name']+'_fold.png'
        fold_hist_fn = os.path.join(infosite_img_path,fold_hist_name)
        fold_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+fold_hist_name
        peak_json[peak_fn]['fold distribution url'] = fold_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)
        hist_stats = [p['fold_enrichment'] for p in hist_peak_stats]
        mp.hist(hist_stats,label=hist_labels,bins=20,log=True)
        mp.title('%s\npeak fold enrichment distribution'%macs_f.file_info['name'])
        mp.xlabel('fold enrichment')
        mp.ylabel('# peaks')
        mp.legend()
        mp.savefig(fold_hist_fn)
        mp.clf()

        fdr_hist_name = macs_f.file_info['name']+'_fdr.png'
        fdr_hist_fn = os.path.join(infosite_img_path,fdr_hist_name)
        fdr_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+fdr_hist_name
        peak_json[peak_fn]['fdr distribution url'] = fdr_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)

        if neg_peak_stats :
            mp.hist(peak_stats['fdr'],label=hist_labels[0],bins=20,log=True)
            mp.title('%s\npeak fdr distribution'%macs_f.file_info['name'])
            mp.xlabel('fdr')
            mp.ylabel('# peaks')
            mp.legend()
            mp.savefig(fdr_hist_fn)
            mp.clf()

        chr_dist_name = macs_f.file_info['name']+'_chr_dist.png'
        chr_dist_fn = os.path.join(infosite_img_path,chr_dist_name)
        chr_dist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+chr_dist_name
        peak_json[peak_fn]['chr distribution url'] = chr_dist_url
        chromos = []
        if json_d.has_key('org') :
            chr_sizes_fn = get_org_settings(json_d['org'])['ucsc_chrom_sizes']
            chromos = [r[0] for r in reader(open(chr_sizes_fn),delimiter='\t')]
        else :
            chromos = set(pos_chr_dist.keys())
            if neg_peak_stats :
                chromos = chromos.union(neg_chr_dist.keys())
            chromos = list(chromos)

        standard_chromos = filter(lambda x: re.search('^chr[0-9MXY]+$',x) is not None,chromos)

        # hack chrM, chrX and chrY so they sort right
        if 'chrM' in standard_chromos :
            standard_chromos[standard_chromos.index('chrM')] = 'chr100'
        if 'chrX' in standard_chromos :
            standard_chromos[standard_chromos.index('chrX')] = 'chr101'
        if 'chrY' in standard_chromos :
            standard_chromos[standard_chromos.index('chrY')] = 'chr102'

        standard_chromos.sort(key=lambda x: int(x.replace('chr','')))

        # unhack chrM, chrX and chrY so they display right
        if 'chr100' in standard_chromos :
            standard_chromos[standard_chromos.index('chr100')] = 'chrM'
        if 'chr101' in standard_chromos :
            standard_chromos[standard_chromos.index('chr101')] = 'chrX'
        if 'chr102' in standard_chromos :
            standard_chromos[standard_chromos.index('chr102')] = 'chrY'

        other_chromos = filter(lambda x: re.search('^chr[0-9MXY]+$',x) is None,chromos)

        pos_plot_chr_dist = defaultdict(int)
        for chrom in standard_chromos :
            pos_plot_chr_dist[chrom] += pos_chr_dist.get(chrom,0)
        for chrom in other_chromos :
            pos_plot_chr_dist['Other'] += pos_chr_dist.get(chrom,0)

        chromos.append('Other')
        mp.figure(figsize=figsize)
        mp.subplots_adjust(bottom=0.18,**subplots_sizes)
        mp.bar(range(len(chromos)),
               [pos_plot_chr_dist[k] for k in chromos],
               width=0.45,
               color='b',
               label='Positive'
              )

        if neg_peak_stats :
            neg_plot_chr_dist = defaultdict(int)
            for chrom in standard_chromos :
                neg_plot_chr_dist[chrom] += neg_chr_dist.get(chrom,0)
            for chrom in other_chromos :
                neg_plot_chr_dist['Other'] += neg_chr_dist.get(chrom,0)

            mp.bar([x+0.45 for x in range(len(chromos))],
                   [neg_plot_chr_dist[k] for k in chromos],
                   width=0.45,
                   color='g',
                   label='Negative'
                  )

        mp.xticks([x+0.45 for x in range(len(chromos))],chromos,rotation=90)
        mp.title('%s\nPeaks by chromosome'%macs_f.file_info['name'])
        mp.xlabel('Chromosome')
        mp.ylabel('# peaks')
        mp.legend()
        mp.savefig(chr_dist_fn)
        mp.clf()

        # pos vs neg peaks
        if neg_peak_stats :
            pos_v_neg_name = '%s_pos_v_neg.png'%macs_f.file_info['name']
            pos_v_neg_fn = os.path.join(infosite_img_path,pos_v_neg_name)
            pos_v_neg_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+pos_v_neg_name
            peak_json[peak_fn]['pos v neg url'] = pos_v_neg_url
            cmd = 'plot_pos_vs_neg_peaks.py --output=%s %s %s'%(pos_v_neg_fn,peak_fn, neg_peak_fn)
            sys.stderr.write(cmd+'\n')
            r = call(cmd,shell=True)

        # motif stuff
        if opts.skip_motif_scan or opts.skip_motif_stuff :
            sys.stderr.write('Obediently skipping motif stuff\n')
        else :
            # not exactly sure the best way to find the filtered macs file yet,
            # just take the .xls file with the longest filename?
            filtered_peak_fns = glob.glob('%s_peaks_*'%macs_f.file_info['name'])
            filtered_peak_fns.sort(key=lambda x: len(x),reverse=True)
            filtered_peak_fn = filtered_peak_fns[0]

            motif_results_fns = glob.glob('%s_motifs_beta*_cv*[0-9].tamo'%macs_f.file_info['name'])
            motif_results_fn = motif_results_fns[0]
            #TODO - do check for file exists

            # motif_scan.py <org> <peak fn> <TAMO motif fn>
            fixed_peak_width = ''
            if json_d['fixed peak width'] != 'none' :
                fixed_peak_width = '--fixed-peak-width=%s'%json_d['fixed peak width']

            cmd = 'motif_scan.py %s --dir=%s/images/ %s %s %s'
            cmd = cmd%(fixed_peak_width,infosite_dir_name,json_d['org'],filtered_peak_fn,motif_results_fn)
            sys.stderr.write(cmd+'\n')
            call(cmd,shell=True)

        # pot_peaks_vs_motifs.py <peaks fn> <seq score fn> <bg score fn>


    # 5. build reSt document
    reSt_fn = exp_name+'_info.rst'
    reSt_path = os.path.join(infosite_path,reSt_fn)
    reSt_html_name = exp_name+'_info.html'
    reSt_html_path = os.path.join(infosite_path,reSt_html_name)
    reSt_url = json_d['stage url'] + '/' + infosite_dir_name + '/' + reSt_html_name
    doc = ReStDocument(reSt_path)
    doc.add(ReStSection("Infopage for %s"%exp_name))

    # basic experiment stats table
    ident = lambda x: x or 'unknown'
    stat_key_labels_fmts = [
                        ('org','Organism',ident),
                        ('analysis path','Analysis Path',ident),
                        ('experiment path','Experiment Path',ident),
                        ('control path','Control Path',ident),
                        ('format','Read Format',ident),
                        ('FDR filter','FDR filter',ident),
                        ('mapping type','Gene Mapping Type',ident),
                        ('mapping window','Gene Mapping Window',lambda x: x and '-%s,%s'%tuple(x)),
                        ('peaks used by THEME','Peaks used by THEME',ident)
                       ]
    stat_rows = [('**%s**'%label, fmt(json_d.get(key))) for key,label,fmt in stat_key_labels_fmts]
    doc.add(ReStSimpleTable(None,stat_rows))

    doc.add(ReStSection('MACS Peak File Stats',level=2))

    # go through peak files
    peak_recs = json_d['peak files']
    fl_str = lambda x: x and '%.2g'%float(x)
    stat_key_labels_fmts = [
                        ('paired peaks','*paired peaks*',ident),
                        ('positive peaks','*positive peaks*',ident),
                        ('negative peaks','*negative peaks*',ident),
                        ('reads under peaks','*reads under positive peaks*',ident),
                        ('total tags in treatment','*Treatment Tags*',ident),
                        ('tags after filtering in treatment','after filtering',ident),
                        ('Redundant rate in treatment','redunancy rate',fl_str),
                        ('maximum duplicate tags at the same position in treatment','max dup. tags',ident),
                        ('total tags in control','*Control Tags*',ident),
                        ('tags after filtering in control','after filtering',ident),
                        ('Redundant rate in control','redunancy rate',fl_str),
                        ('maximum duplicate tags at the same position in control','max dup. tags',ident),
                        ('peak tag count filter','*Minimum peak tag count*',ident),
                        ('d','*MACS d*',ident),
                        ('band width','*band width*',ident),
                        ('MACS version','*MACS version*',ident),
                        ('pvalue cutoff','*p-value cutoff*',lambda x: '1e%d'%int(log(x,10))),
                       ]

    for peak_fn,peak_stats in peak_recs.items() :

        # add the new section and stats table
        doc.add(ReStSection(peak_fn,level=3))
        stat_rows = [('*%s*'%label, fmt(peak_stats.get(key))) for key,label,fmt in stat_key_labels_fmts]
        doc.add(ReStSimpleTable(None,stat_rows))

        # link to the peaks file
        peak_infosite_name = os.path.join(infosite_dir_name,peak_fn)
        peak_infosite_path = os.path.abspath(peak_infosite_name)
        peak_infosite_url = json_d['stage url'] + '/' + peak_infosite_name
        call('cp %s %s'%(peak_fn,os.path.join(infosite_dir_name,peak_fn)),shell=True)
        doc.add(ReStSimpleTable(None,[('**MACS Peaks File**','`%s`_'%peak_infosite_url)]))
        doc.add(ReStHyperlink(peak_infosite_url,url=peak_infosite_url))

        # UCSC track info
        if peak_stats.has_key('ucsc tracks') :
            ucsc_tbl = ReStSimpleTable(('**UCSC Genome Browser Track Lines**',),
                                      [[x] for x in peak_stats['ucsc tracks']])
            doc.add(ucsc_tbl)
        else :
            doc.add(ReStSimpleTable(None,[['UCSC integration was not enabled for this experiment']]))

        # peak quality plots
        if neg_peak_stats :
            img_tbl1 = ReStSimpleTable(None, [
                        [
                         ReStImage(peak_stats['pos v neg url'],options={'width':'600px','align':'center'}),
                        ]
                       ]
                      )
            doc.add(img_tbl1)

            fdr_val = ReStImage(peak_stats['fdr distribution url'],options={'width':'250px','align':'center'})
        else :
            fdr_val = 'No negative peaks, not applicable'

        img_tbl2 = ReStSimpleTable(None, [
                    [
                     ReStImage(peak_stats['length distribution url'],options={'width':'250px','align':'center'}),
                     ReStImage(peak_stats['tag distribution url'],options={'width':'250px','align':'center'}),
                     ReStImage(peak_stats['pvalue distribution url'],options={'width':'250px','align':'center'})
                    ],
                    [
                     ReStImage(peak_stats['fold distribution url'],options={'width':'250px','align':'center'}),
                     fdr_val,
                     ReStImage(peak_stats['chr distribution url'],options={'width':'250px','align':'center'})
                    ]
                  ]
                  )
        doc.add(img_tbl2)

        # gene info
        gene_fn = peak_stats['name']+'_genes.txt'
        gene_link = os.path.join(infosite_dir_name,gene_fn)
        if not os.path.exists(gene_link) :
            shutil.copyfile(gene_fn,gene_link)
        gene_url = json_d['stage url']+'/'+gene_link

        # gather other gene mapping stats
        # knownGeneID
        # geneSymbol
        # chr
        # start
        # end
        # length
        # summit
        # tags
        # -10*log10(pvalue)
        # fold_enrichment
        # FDR(%)
        # peak
        # loc
        # dist
        # from
        # feature
        # score
        # map
        # type
        # map
        # subtype

        gene_reader = DictReader(open(gene_fn),delimiter='\t')
        gene_stats = defaultdict(set)
        gene_pvals = defaultdict(float)
        for rec in gene_reader :
            gene_stats['num knownGenes'].add(rec['knownGeneID'])
            gene_stats['num geneSymbols'].add(rec['geneSymbol'])
            gene_pvals[rec['geneSymbol']] = max(gene_pvals[rec['geneSymbol']],float(rec['-10*log10(pvalue)']))
        gene_pvals = gene_pvals.items()
        gene_pvals.sort(key=lambda x: x[1],reverse=True)
        for k,v in gene_pvals[:20]:
            print k,v
        gene_mapping_data = [('**# knownGenes mapped**',len(gene_stats['num knownGenes'])),
                             ('**# gene symbols mapped**',len(gene_stats['num geneSymbols'])),
                             ('**Top 10 gene symbols**',','.join([x[0] for x in gene_pvals[:10]])),
                             ('**All gene mappings**','`%s`_'%gene_url)
                            ]

        # plots from plot_peak_loc_dist.py
        gene_pie_name = exp_name+'_gene_map.png'
        peak_pie_name = exp_name+'_peak_map.png'
        hist_name = exp_name+'_peak_dist.png'
        pval_bar_name = exp_name+'_pval_bar.png'
        peak_loc_d = {'out_dir':infosite_path,
                      'gene_pie_fn':os.path.join(infosite_path,'images',gene_pie_name),
                      'peak_pie_fn':os.path.join(infosite_path,'images',peak_pie_name),
                      'pval_bar_fn':os.path.join(infosite_path,'images',pval_bar_name),
                      'hist_fn':os.path.join(infosite_path,'images',hist_name),
                      'peak_fn':peak_fn,
                      'gene_name':gene_fn
                      }
        cmd = 'plot_peak_loc_dist.py --save -d %(out_dir)s -g %(gene_pie_fn)s ' \
              '-p %(peak_pie_fn)s -f %(hist_fn)s -b %(pval_bar_fn)s ' \
              '%(peak_fn)s %(gene_name)s'
        sys.stderr.write(cmd%peak_loc_d+'\n')
        call(cmd%peak_loc_d,shell=True)
        peak_stats['gene map url'] = json_d['stage url']+'/'+infosite_dir_name+'/images/'+gene_pie_name
        peak_stats['peak map url'] = json_d['stage url']+'/'+infosite_dir_name+'/images/'+peak_pie_name
        peak_stats['pval bar url'] = json_d['stage url']+'/'+infosite_dir_name+'/images/'+pval_bar_name
        peak_stats['dist url'] = json_d['stage url']+'/'+infosite_dir_name+'/images/'+hist_name

        # make links to the different peaks files
        feature_patts = ('promoter.txt','gene_exon.txt','gene_intron.txt','after.txt','intergenic.xls')
        feature_data = []
        feature_urls = []

        for patt in feature_patts :
            feature_fn = '%s_*_%s'%(peak_stats['name'],patt)
            feature_path = glob.glob(os.path.join(infosite_dir_name,feature_fn))
            if len(feature_path) == 0 :
                sys.stderr.write('Warning: %s could not be found, skipping feature type\n'%os.path.join(infosite_dir_name,feature_fn))
                continue
            feature_path = feature_path[0]
            feature_url = json_d['stage url']+'/'+feature_path

            # create UCSC formatted versions of the files
            if patt.endswith('.txt') : # these have gene columns
                feature_type = patt.replace('.txt','')
                ucsc_feature_fn = feature_fn.replace('.txt','_ucsc.txt')
                st,en = 2,4
            elif patt.endswith('.xls') :
                feature_type = patt.replace('.xls','')
                ucsc_feature_fn = feature_fn.replace('.xls','_ucsc.xls')
                st,en = 0,2

            ucsc_feature_path = os.path.join(infosite_dir_name,ucsc_feature_fn)
            ucsc_feature_f = open(ucsc_feature_path,'w')
            ucsc_feature_writer = writer(ucsc_feature_f,delimiter='\t')
            for l in reader(open(feature_path),delimiter='\t') :
                rec = l[0:st] + \
                      ['%s:%s-%s'%tuple(l[st:en+1])] + \
                      l[en+1:]
                ucsc_feature_writer.writerow(rec)
            ucsc_feature_f.close()

            ucsc_feature_url = json_d['stage url']+'/'+ucsc_feature_path

            feature_data.append(('**%s peaks**'%feature_type,'`%s`_ `UCSC %s`_'%(feature_url,feature_type)))
            feature_urls.append(ReStHyperlink(feature_url,url=feature_url))
            feature_urls.append(ReStHyperlink('UCSC %s'%feature_type,url=ucsc_feature_url))

        gene_mapping_data.extend(feature_data)
        feat_tbl = ReStSimpleTable(('**Gene mapping data**',''),gene_mapping_data)
        doc.add(feat_tbl)
        doc.add(ReStHyperlink(gene_url,url=gene_url))
        for url in feature_urls :
            doc.add(url)

        img_tbl3 = ReStSimpleTable(None, [
                    [
                     ReStImage(peak_stats['gene map url'],options={'align':'center'}),
                     ReStImage(peak_stats['peak map url'],options={'align':'center'})
                    ],
                    [
                     ReStImage(peak_stats['pval bar url'],options={'align':'center'}),
                     ReStImage(peak_stats['dist url'],options={'align':'center'})
                    ]
                   ]
                   )
        doc.add(img_tbl3)

        # now put some motif stuff up there


        if opts.skip_motif_stuff :
            sys.stderr.write('Obediently skipping even more motif stuff\n')
        else :
            # THEME refines all motifs, display the top 30

            # for now, just list a table of the top 30 significant, unrefined motifs
            doc.add(ReStSection('%s Top 30 Refined Motif Results'%peak_stats['name'],level=3))
            motif_results_fns = glob.glob('%s_motifs_beta*_cv*[0-9].txt'%macs_f.file_info['name']) #catRun_mfold10,30_pval1e-5_motifs_beta0.0_cv5.txt
            #TODO - do check for file exists

            motif_results_fn = motif_results_fns[0]

            motif_reader = reader(open(motif_results_fn),delimiter='\t')

            motif_header = motif_reader.next()
            motif_data = []
            top_n = 30
            motif_fmts = (ident,ident,int,fl_str,fl_str,fl_str,fl_str,fl_str,fl_str)
            motif_plot_urls = []
            for rec in motif_reader :
                motif_data.append([f(x) for f,x in zip(motif_fmts,rec)])
                """
                if rec[2] in motif_sig_inds_d.keys() :
                    from_id = motif_sig_inds_d[rec[2]]
                    try :
                        old_id_fn = glob.glob(infosite_dir_name+'/images/*_%d_peakmot.png'%from_id)[0]
                        new_id_fn = old_id_fn.replace('_%d_'%from_id,'_%s_'%rec[2])
                        os.rename(old_id_fn,new_id_fn)
                    except :
                        sys.stderr.write("Couldn't rename file for pattern %s, just " \
                                         "assuming its there\n"%(infosite_dir_name+'/images/*_%d_peakmot.png'%from_id))
                """
                try :
                    new_id_fn = glob.glob(infosite_dir_name+'/images/*_%s_peakmot.png'%rec[2])[0]
                    motif_plot_urls.append(json_d['stage url']+'/'+new_id_fn)
                except Exception, e :
                    sys.stderr.write('Exception thrown in plotting motifs, skipping: %s'%e)

            doc.add(ReStSimpleTable(['**%s**'%x for x in motif_header],motif_data[:top_n]))

            # create another file with the full table
            motif_results_base, motif_results_ext = os.path.splitext(motif_results_fn)
            motif_doc_fn = motif_results_base+'.rst'
            motif_doc_path = os.path.join(infosite_path,motif_doc_fn)
            motif_doc_html_fn = motif_results_base+'.html'
            motif_doc_html_path = os.path.join(infosite_path,motif_doc_html_fn)
            motif_doc_url = json_d['stage url']+'/'+infosite_dir_name+'/'+motif_doc_html_fn
            motif_doc = ReStDocument(motif_doc_path)
            motif_doc.add(ReStSection('%s Full Motif Results'%peak_stats['name']))
            motif_doc.add('`Back to main infopage`_')
            motif_doc.add(ReStSimpleTable(['**%s**'%x for x in motif_header],motif_data))
            motif_doc.add('`Back to main infopage`_')
            motif_doc.add(ReStHyperlink('Back to main infopage',url=reSt_url))
            motif_doc.write()
            motif_doc.close()
            rst2html_call = 'rst2html.py --stylesheet-path=/nfs/antdata/web_stage/css/lsr.css ' \
                            '%s %s'%(motif_doc_path,motif_doc_html_path)
            sys.stderr.write(rst2html_call+'\n')
            r = call(rst2html_call,shell=True)
            doc.add('`All refined motifs`_')
            doc.add(ReStHyperlink('All refined motifs',url=motif_doc_url))

            # individual motif plots
            plt_tbl = []
            for i,url in enumerate(motif_plot_urls[:30]) :
                if i%3 == 0 :
                    plt_tbl.append([])
                plt_tbl[-1].append(ReStImage(url))

            doc.add(ReStSimpleTable(('**Peak strength vs refined motif strength**','(based on top 2000 peak sequences by pvalue)',''),plt_tbl))

    doc.write()
    doc.close()

    # 6. convert reSt to PDF and HTML
    rst2html_call = 'rst2html.py --stylesheet-path=/nfs/antdata/web_stage/css/lsr.css ' \
                    '%s %s'%(reSt_path,reSt_html_path)
    sys.stderr.write(rst2html_call+'\n')
    r = call(rst2html_call,shell=True)

    pdf_name = exp_name+'_info.pdf'
    pdf_path = os.path.join(infosite_path,pdf_name)
    r = call('rst2pdf %s -o %s'%(reSt_path,pdf_path),shell=True)

    # 7. write out url to infosite
    print json_d['stage url']+'/'+infosite_dir_name+'/'+reSt_html_name
    open(infosite_dir_name+'_url.txt','w').write(json_d['stage url']+'/'+infosite_dir_name+'/'+reSt_html_name+'\n')
