#!/usr/bin/env python

import glob
import json
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as mp
import os
import sys

from collections import defaultdict
from math import log
from optparse import OptionParser
from subprocess import call

from chipsequtil import MACSFile
from reStUtil import *

usage = '%prog [options] [<peak filename> <peak filename> ...]'
parser = OptionParser(usage=usage)
parser.add_option('-d','--dir',dest='dir',default='.',help='Source directory [default: %default]')
parser.add_option('-n','--name',dest='name',help='Experiment name [default: current directory name]')

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
                      'stage path':'/nfs/antdata/web_stage/'+curr_user
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
        macs_f = MACSFile(peak_fn)
        peak_json[peak_fn] = macs_f.file_info

        peak_stats = defaultdict(list)
        num_peaks = 0
        for peak in macs_f :
            peak_stats['length'].append(peak['length'])
            peak_stats['tags'].append(peak['tags'])
            peak_stats['pvalue'].append(peak['-10*log10(pvalue)'])
            peak_stats['fold_enrichment'].append(peak['fold_enrichment'])
            peak_stats['fdr'].append(peak['FDR(%)'])
            num_peaks += 1

        peak_json[peak_fn]['number of peaks'] = num_peaks

        font = {'size':'9'}
        mp.rc('font',**font)

        figsize = (3.5,3.5)
        subplots_sizes = {'top':0.8,'left':0.15,'right':0.95}
        # create histograms for each of the attributes
        len_hist_name = macs_f.file_info['name']+'_length.png'
        len_hist_fn = os.path.join(infosite_img_path,len_hist_name)
        len_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+len_hist_name
        peak_json[peak_fn]['length distribution url'] = len_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)
        mp.hist(peak_stats['length'],bins=20,log=True)
        mp.title('%s\npeak length distribution'%macs_f.file_info['name'])
        mp.xlabel('peak length')
        mp.ylabel('# peaks')
        mp.savefig(len_hist_fn)
        mp.clf()

        tags_hist_name = macs_f.file_info['name']+'_tags.png'
        tags_hist_fn = os.path.join(infosite_img_path,tags_hist_name)
        tags_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+tags_hist_name
        peak_json[peak_fn]['tag distribution url'] = tags_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)
        mp.hist(peak_stats['tags'],bins=20,log=True)
        mp.title('%s\npeak tag count distribution'%macs_f.file_info['name'])
        mp.xlabel('# tags')
        mp.ylabel('# peaks')
        mp.savefig(tags_hist_fn)
        mp.clf()

        pval_hist_name = macs_f.file_info['name']+'_pval.png'
        pval_hist_fn = os.path.join(infosite_img_path,pval_hist_name)
        pval_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+pval_hist_name
        peak_json[peak_fn]['pvalue distribution url'] = pval_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)
        mp.hist(peak_stats['pvalue'],bins=20,log=True)
        mp.title('%s\npeak -10*log10(p-valuek) distribution'%macs_f.file_info['name'])
        mp.xlabel('-10*log10(p-value)')
        mp.ylabel('# peaks')
        mp.savefig(pval_hist_fn)
        mp.clf()

        fold_hist_name = macs_f.file_info['name']+'_fold.png'
        fold_hist_fn = os.path.join(infosite_img_path,fold_hist_name)
        fold_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+fold_hist_name
        peak_json[peak_fn]['fold distribution url'] = fold_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)
        mp.hist(peak_stats['fold_enrichment'],bins=20,log=True)
        mp.title('%s\npeak fold enrichment distribution'%macs_f.file_info['name'])
        mp.xlabel('fold enrichment')
        mp.ylabel('# peaks')
        mp.savefig(fold_hist_fn)
        mp.clf()

        fdr_hist_name = macs_f.file_info['name']+'_fdr.png'
        fdr_hist_fn = os.path.join(infosite_img_path,fdr_hist_name)
        fdr_hist_url = json_d['stage url']+'/'+infosite_dir_name+'/images/'+fdr_hist_name
        peak_json[peak_fn]['fdr distribution url'] = fdr_hist_url
        mp.figure(figsize=figsize)
        mp.subplots_adjust(**subplots_sizes)
        mp.hist(peak_stats['fdr'],bins=20,log=True)
        mp.title('%s\npeak fdr distribution'%macs_f.file_info['name'])
        mp.xlabel('fdr')
        mp.ylabel('# peaks')
        mp.savefig(fdr_hist_fn)
        mp.clf()

    # 5. build reSt document
    reSt_fn = exp_name+'_info.rst'
    reSt_path = os.path.join(infosite_path,reSt_fn)
    doc = ReStDocument(reSt_path)
    doc.add(ReStSection("Infopage for %s"%exp_name))

    # basic experiment stats table
    stat_header = ('','')
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
    doc.add(ReStSimpleTable(stat_header,stat_rows))

    doc.add(ReStSection('MACS Peak File Stats',level=2))

    # go through peak files
    peak_fns = json_d['peak files']
    fl_str = lambda x: x and '%.2f'%float(x)
    stat_key_labels_fmts = [
                        ('total tags in treatment','*Treatment Tags*',ident),
                        ('tags after filtering in treatment','after filtering',ident),
                        ('Redundant rate in treatment','redunancy rate',fl_str),
                        ('maximum duplicate tags at the same position in treatment','max dup. tags',ident),
                        ('total tags in control','*Control Tags*',ident),
                        ('tags after filtering in control','after filtering',ident),
                        ('Redundant rate in control','redunancy rate',fl_str),
                        ('maximum duplicate tags at the same position in control','max dup. tags',ident),
                        ('d','*MACS d*',ident),
                        ('band width','*band width*',ident),
                        ('MACS version','*MACS version*',ident),
                        ('pvalue cutoff','*p-value cutoff*',lambda x: '1e%d'%int(log(x,10))),
                        ('number of peaks','*number of peaks*',ident),
                        #('experiment path','Experiment Path',ident),
                        #('control path','Control Path',ident),
                        #('format','Read Format',ident),
                        #('FDR filter','FDR filter',ident),
                        #('mapping type','Gene Mapping Type',ident),
                        #('mapping window','Gene Mapping Window',lambda x: '-%s,%s'%tuple(x)),
                        #('peaks used by THEME','Peaks used by THEME',ident)
                       ]

    for peak_fn,peak_stats in peak_fns.items() :
        doc.add(ReStSection(peak_fn,level=3))
        stat_rows = [('*%s*'%label, fmt(peak_stats.get(key))) for key,label,fmt in stat_key_labels_fmts]
        doc.add(ReStSimpleTable(('**Peak stats**',''),stat_rows))
        doc.add(ReStImage(peak_stats['length distribution url']))
        doc.add(ReStImage(peak_stats['tag distribution url']))
        doc.add(ReStImage(peak_stats['pvalue distribution url']))
        doc.add(ReStImage(peak_stats['fold distribution url']))
        doc.add(ReStImage(peak_stats['fdr distribution url']))

    # gene info

    # now put some motif stuff up there

    doc.write()
    doc.close()

    # 6. convert reSt to PDF and HTML
    html_name = exp_name+'_info.html'
    html_path = os.path.join(infosite_path,html_name)
    rst2html_call = 'rst2html --stylesheet-path=/nfs/antdata/web_stage/css/lsr.css ' \
                    '%s %s'%(reSt_path,html_path)
    r = call(rst2html_call,shell=True)

    pdf_name = exp_name+'_info.pdf'
    pdf_path = os.path.join(infosite_path,pdf_name)
    #r = call('rst2pdf %s %s'%(reSt_path,pdf_path),shell=True)

    # 7. write out url to infosite
    print json_d['stage url'] + '/' + infosite_dir_name + '/' + html_name
