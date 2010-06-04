import math
import os

from ConfigParser import ConfigParser
from csv import DictReader
from collections import defaultdict

import chipsequtil

# for RefGeneDB
from util import KeyedBinaryTree


def get_file_parts(path) :
    """For <path>/<basename>.<ext>, returns 4-tuple (<path>,<basename>.<ext>,<basename>,<ext>)"""
    path,fn = os.path.split(path)
    basename,ext = os.path.splitext(fn)
    return path,fn,basename,ext

def parse_number(n) :
    """Try to cast intput first to float, then int, returning unchanged if both fail"""
    try :
        return float(n) if '.' in n else int(n)
    except :
        return n

class GERALDOutput :
    """Container for one line of GERALD alignment output"""

    FIELD_NAMES = ['machine',
                   'run_number',
                   'lane',
                   'tile',
                   'x_coord',
                   'y_coord',
                   'index',
                   'read_no',
                   'read',
                   'quality_string',
                   'match_chromo',
                   'match_contig',
                   'match_pos',
                   'match_strand',
                   'match_desc',
                   'single_read_score',
                   'paired_read_score',
                   'partner_chromo',
                   'partner_contig',
                   'partner_offset',
                   'partner_strand',
                   'filtering',
                   ]

    def __init__(self,line) :

        if type(line) == str :
            line = line.strip().split('\t')

        if len(line) != len(GERALDOutput.FIELD_NAMES) :
            raise GERALDOutput.FormatException('Expected %d fields in input, \
                                               found %d in line: %s'%
                                               (len(GERALDOutput.FIELD_NAMES),
                                                len(line),
                                                line))

        for fn,d in zip(GERALDOutput.FIELD_NAMES,line) :
            setattr(self,fn,parse_number(d))

    def __repr__(self) :
        return 'GERALDOutput(%s)'%repr(self.output_format())

    def output_format(self) :
        """Tab delimited string of fields as they would appear in GERALD output file"""
        return '\t'.join([str(getattr(self,d)) for d in GERALDOutput.FIELD_NAMES])+'\n'

    class FormatException(Exception) :
        """GERALD format exception, raised on malformatted input"""
        pass

class BEDOutput :
    """Container for one line of BED alignment output"""

    FIELD_NAMES = ['chrom',
                   'chromStart',
                   'chromEnd',
                   'name',
                   'score',
                   'strand',
                   'thickStart',
                   'thickEnd',
                   'itemRgb',
                   'blockCount',
                   'blockSizes',
                   'blockStarts',
                   ]

    def __init__(self,line='',*args,**kwargs) :

        if type(line) == str :
            line = line.strip().split('\t')

        if len(line) < 3 and any([x not in kwargs.keys() for x in ['chrom','chromStart','chromEnd']]) :
            raise BEDOutput.FormatException('Format requres at least 3 fields in \
                                            input, found %d in line: %s'%(len(line),line))
        if len(line) > len(BEDOutput.FIELD_NAMES) :
            raise BEDOutput.FormatException('Format requres at most %d fields in \
                                             input, found %d in line: %s'%
                                             (len(BEDOutput.FIELD_NAMES),len(line),line))

        empty_fields = ['']*(len(BEDOutput.FIELD_NAMES)-len(line))
        for fn,d in zip(BEDOutput.FIELD_NAMES,line+empty_fields) :
            setattr(self,fn,parse_number(d))

        # kwargs override line input
        for k,v in kwargs.items() :
            setattr(self,k,parse_number(v))

    def __repr__(self) :
        return 'BEDOutput(%s)'%(repr(self.output_format()))

    def output_format(self) :
        """Returns a string for the BED line as it would appear in a file"""
        return '\t'.join([str(getattr(self,d)) for d in BEDOutput.FIELD_NAMES])+'\n'

    class FormatException(Exception) :
        """BED format exception, raised on malformatted input"""
        pass


class BEDFile(DictReader) :
    '''An iterable object (subclasses csv.DictReader) containing the records in
    the supplied BED formatted file'''
    def __init__(self,bed_fn) :
        DictReader.__init__(self,open(bed_fn),delimiter='\t',
                            fieldnames=BEDOutput.FIELD_NAMES)


class RefGeneOutput(object) :
    # http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.sql
    FIELD_NAMES = ['bin',
                   'name',
                   'chrom',
                   'strand',
                   'txStart',
                   'txEnd',
                   'cdsStart',
                   'cdsEnd',
                   'exonCount',
                   'exonStarts',
                   'exonEnds',
                   'score',
                   'name2',
                   'cdsStartStat',
                   'cdsEndStat',
                   'exonFrames',]


class RefGeneFile(DictReader) :
    '''An iterable object (subclasses csv.DictReader) containing the records in
    the supplied BED formatted file'''
    def __init__(self,refGene_fn) :
        refGene_f = open(refGene_fn)
        # check for header
        first_line = refGene_f.next()
        if not first_line.strip().startswith('#') :
            refGene_f.seek(0) # first line not header, reset the file pointer
        DictReader.__init__(self,refGene_f,delimiter='\t',fieldnames=RefGeneOutput.FIELD_NAMES)


#TODO maybe, finish this
class RefGeneDB :
    '''A class for querying RefGene annotation files. NOT DONE.'''

    def __init__(self,refgene_fn) :
        self._chrom_trees = defaultdict(KeyedBinaryTree)
        refgene_f = RefGeneFile(refgene_fn)
        genes = defaultdict(list)
        for gene in refgene_f :
            genes[gene['chrom']].append(gene)

        # do stuff to ensure a balanced tree for each chromosome
        for chrom,gene_list in genes.items() :
            gene_list.sort(key=lambda x: int(x['txStart']))
            first_half, second_half = gene_list[:len(gene_list)/2],gene_list[len(gene_list)/2:]
            first_half.reverse()
            for i in range(min(len(first_half,second_half))) :
                to_add = first_half.pop(i)
                self._chrom_trees[chrom].addNode(int(to_add['txStart']),to_add)

def gerald_to_bed(gerald,min_fields=False) :
    """Convert a GERALDOutput object into a BEDOutput object

    Keyword argument *min_fields* produces BED alignment with only the first 
    three fields populated
    """

    d = {}.fromkeys(BEDOutput.FIELD_NAMES,'')

    # required BED fields
    d['chrom'] = gerald.match_chromo
    d['chromStart'] = gerald.match_pos
    d['chromEnd'] = gerald.match_pos+len(gerald.read)

    # load the remaining information
    if not min_fields :
        d['strand'] = '+' if gerald.match_strand == 'F' else '-'
        # TODO consider encoding single-read alignment score into BED score format
        # that's it?
    return BEDOutput(**d)


class MACSOutput(object) :
    FIELD_NAMES = ['chr',
                   'start',
                   'end',
                   'length',
                   'summit',
                   'tags',
                   '-10*log10(pvalue)',
                   'fold_enrichment',
                   'FDR(%)',
                  ]

class MACSFile(DictReader) :
    '''An iterable object (subclasses csv.DictReader) containing the records in
    the supplied MACS formatted peak file'''
    def __init__(self,macs_fn) :
        self.meta_data = []
        self.file_info = {}
        f = open(macs_fn)
        done_with_header = False
        while not done_with_header :
            l = f.next().strip()
            if l.startswith('#') :
                if l.count('=') != 0 :
                    k,v = l.split('=',1)
                    self.file_info[k.strip()] = parse_number(v.strip())
                self.meta_data.append(l)
            elif l == '\t'.join(MACSOutput.FIELD_NAMES) :
                self.meta_data.append(l)
                done_with_header = True

        DictReader.__init__(self,f,delimiter='\t',fieldnames=MACSOutput.FIELD_NAMES)



GLOBAL_SETTINGS_FN = os.path.join(os.path.split(chipsequtil.__file__)[0],'org_settings.cfg')
LOCAL_SETTINGS_FN = os.path.expanduser(os.path.join('~','.org_settings.cfg'))
_ALL_SETTINGS, _LOCAL_SETTINGS, _GLOBAL_SETTINGS = range(3)

def _get_org_settings(org_key=None,addnl_configs=[],src=_ALL_SETTINGS) :
    """Utility function used by get_org_settings and get_all_settings, should \
    not be called directly"""

    config = ConfigParser()
    chipsequtil_base =     conf_fns = []
    if src in [_LOCAL_SETTINGS, _ALL_SETTINGS] :
        conf_fns.append(LOCAL_SETTINGS_FN)
    if src in [_GLOBAL_SETTINGS, _ALL_SETTINGS] :
        conf_fns.append(GLOBAL_SETTINGS_FN)
    config.read(conf_fns+addnl_configs)

    d = {}
    if org_key is None :
        for sec in config.sections() :
            # try to cast numeric-looking arguments into float, int
            d[sec] = dict([(k,parse_number(v)) for k,v in config.items(sec)])
    else :
        d = dict([(k,parse_number(v)) for k,v in config.items(org_key)])

    return d


def get_org_settings(org_key,addnl_configs=[]) :
    '''Returns a dict of setting/path values for a given organism as specified in system-wide and user's settings.'''
    return _get_org_settings(org_key,addnl_configs=[])


def get_all_settings(addnl_configs=[]) :
    '''Returns a dict of setting/path values for every organism as specified in system-wide and user's settings.'''
    return _get_org_settings(None,addnl_configs=addnl_configs)


def get_global_settings() :
    '''Returns a dict of the global setting/path values'''
    return _get_org_settings(None,src=_GLOBAL_SETTINGS)


def get_local_settings() :
    '''Returns a dict of the current user's setting/path values'''
    return _get_org_settings(None,src=_LOCAL_SETTINGS)


RC_MAP_TABLE = ''.join([chr(i) for i in xrange(256)])
# use pairs to replace, otherwise we overwrite previous substitutions
RC_MAP_TABLE = RC_MAP_TABLE.replace('ab','tb').replace('tu','au').replace('cd','gd').replace('gh','ch')
RC_MAP_TABLE = RC_MAP_TABLE.replace('AB','TB').replace('TU','AU').replace('CD','GD').replace('GH','CH')
def reverse_complement(seq) :
    """Reverse complements nucleotide string *seq*.  Leaves non-nucleotide characters uneffected."""
    seq_complement = list(seq.translate(RC_MAP_TABLE))
    seq_complement.reverse()
    return ''.join(seq_complement)


def get_gc_content(seq) :
    '''returns the GC content of a DNA sequence as python string'''
    seq = seq.lower()
    return (seq.count('c')+seq.count('g'))/float(len(seq))


def get_gc_content_distribution(sequences,bins=100) :
    '''returns a function that approximates the GC content distibution of the
    provided sequences.  Approximation is performed by binning.'''
    gc_contents = [get_gc_content(s) for s in sequences]
    gc_contents.sort()

    # count up the sequences for each bin
    bin_counts = [0.]*bins
    for c in gc_contents :
        sample_bin = int(math.floor(c*bins))
        bin_counts[sample_bin] += 1

    # normalize bin counts
    norm_bins = [x/len(sequences) for x in bin_counts]

    # create a closure for this set of sequences
    #def f(seq) :
    #    gc = get_gc_content(seq)
    #    return norm_bins[int(math.floor(gc*bins))]

    return norm_bins


def get_size_distribution(sequences) :
    return (len(s) for s in sequences)



