import math
import os
import re
import sys

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


class SmartFileIter :

    def __init__(self,f) :
        if not hasattr(self,'FIELD_NAMES') or not hasattr(self,'FIELD_TYPES') :
            raise Exception('Subclasses must define class members FIELD_NAMES and FIELD_TYPES')
        if isinstance(f,str) :
            f = open(f)
        self._dict_reader = DictReader(f,delimiter='\t',fieldnames=self.FIELD_NAMES)
        self.fieldnames = self.FIELD_NAMES
        self.curr_line = self._dict_reader.next()

        if all([self.curr_line[k] is None for k in self.FIELD_NAMES[1:]]) : # probably space delimited, try again
            f.seek(0)
            self._dict_reader = DictReader(f,delimiter=' ',fieldnames=self.FIELD_NAMES)
            self.curr_line = self._dict_reader.next()

        if self.FIELD_NAMES[0] in self.curr_line.values() :
            self.curr_line = self._dict_reader.next()

    def __iter__(self) :
        return self

    def __getattr__(self,attr) :
        try:
            return self.__dict__[attr]
        except KeyError :
            return getattr(self._dict_reader,attr)

    def next(self) :

        if self.curr_line is None :
            raise StopIteration()

        line = self.curr_line
        for k,f in zip(self.FIELD_NAMES, self.FIELD_TYPES) :
            try :
                line[k] = f(line[k])
            except Exception, e :
                #sys.stderr.write('Warning: field %s on line %d could not be properly formatted, exception %s\n'%(k,self._dict_reader.reader.line_num,str(e)))
                line[k] = line[k]

        try :
            self.curr_line = self._dict_reader.next()
        except StopIteration :
            self.curr_line = None

        return line


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


class BEDFile(SmartFileIter) :
    FIELD_NAMES = BEDOutput.FIELD_NAMES
    FIELD_TYPES = [str,int,int,str,float,str,int,int,str,lambda x: x.split(','), lambda x: x.split(','), lambda x: x.split(',')]


class BEDFile_dictreader(DictReader) :
    '''An iterable object (subclasses csv.DictReader) containing the records in
    the supplied BED formatted file.'''
    FIELD_NAMES = BEDOutput.FIELD_NAMES
    def __init__(self,bed) :
        '''*bed* is either a filename or a file-like object representing a BED file'''
        if isinstance(bed,str) :
            bed = open(bed)
        DictReader.__init__(self,bed,delimiter='\t',
                            fieldnames=BEDOutput.FIELD_NAMES)


class AffyBiocFile(DictReader) :
    '''An iterable object (subclasses csv.DictReader) containing microarray data records in
    the supplied bioconductor formatted file.'''

    FIELD_NAMES = [ 'ID',
                    'Symbol',
                    'Name',
                    'M',
                    'A',
                    't',
                    'P.Value',
                    'B'
                  ]

    def __init__(self,affyfn) :
        '''*affyfn* is either a filename or a file-like object representing a bioconductor output file'''
        if isinstance(affyfn,str) :
            bed = open(bed)
        DictReader.__init__(self,bed,delimiter='\t',
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

class KnownGeneFile :

    FIELD_NAMES = [ 'name',
                    'chrom',
                    'strand',
                    'txStart',
                    'txEnd',
                    'cdsStart',
                    'cdsEnd',
                    'exonCount',
                    'exonStarts',
                    'exonEnds',
                    'proteinID',
                    'alignID',
                  ]

    # function pointers for correct formatting of field names
    FIELD_TYPES = [ str,
                    str,
                    str,
                    int,
                    int,
                    int,
                    int,
                    lambda x: [int(y) for y in x.split(',') if len(y) > 0],
                    lambda x: [int(y) for y in x.split(',') if len(y) > 0],
                    lambda x: [int(y) for y in x.split(',') if len(y) > 0],
                    str,
                    str,
                  ]

    def __init__(self,kg_fn) :
        self.meta_data = []
        self.file_info = {}
        f = open(kg_fn)
        self._dict_reader = DictReader(f,delimiter='\t',fieldnames=KnownGeneFile.FIELD_NAMES)

    def __iter__(self) :
        return self

    def next(self) :
        line = self._dict_reader.next()
        for k,f in zip(self.FIELD_NAMES,self.FIELD_TYPES) :
            line[k] = f(line[k])
        return line


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


class MACSFile(SmartFileIter) :
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

    FIELD_TYPES = [str,
                   int,
                   int,
                   int,
                   int,
                   int,
                   float,
                   float,
                   float
                  ]

    _METADATA_REGEXES = [
            u'# This file is generated by (MACS version) (.*)',
            u'# (name) = (.*)',
            u'# (format) = (.*)',
            u'# (ChIP-seq file) = (.*)',
            u'# (control file) = (.*)',
            u'# (effective genome size) = (.*)',
            u'# (band width) = (\d+)',
            u'# (model fold) = (.*)',
            u'# (pvalue cutoff) = (.*)',
            u'# (Range for calculating regional lambda) is: (.*)',
            u'# (tag size) is determined as (\d+) bps',
            u'# (total tags in treatment): (\d+)',
            u'# (tags after filtering in treatment): (\d+)',
            u'# (maximum duplicate tags at the same position in treatment) = (\d+)',
            u'# (Redundant rate in treatment): (.*)',
            u'# (total tags in control): (.*)',
            u'# (tags after filtering in control): (.*)',
            u'# (maximum duplicate tags at the same position in control) = (\d+)',
            u'# (Redundant rate in control): (.*)',
            u'# (d) = (\d+)'
            ]

    '''An iterable object (subclasses csv.DictReader) containing the records in
    the supplied MACS formatted peak file'''
    def __init__(self,macs_f) :
        self.meta_data = []
        self.file_info = {}
        if isinstance(macs_f,str) :
            f = open(macs_f)
        else :
            f = macs_f
        done_with_header = False
        while not done_with_header :
            l = f.next().strip()
            if l.startswith('#') :
                for regex in MACSFile._METADATA_REGEXES :
                    m = re.search(regex,l)
                    if m is not None :
                        self.file_info[m.group(1).strip()] = parse_number(m.group(2).strip())
                self.meta_data.append(l)
            elif l.startswith('\t'.join(MACSOutput.FIELD_NAMES[:5])) :
                self.meta_data.append(l)
                done_with_header = True

        SmartFileIter.__init__(self,f)


# for backwards compatibility, use MACSFile instead...?
class MACSOutput(object) :
    FIELD_NAMES = MACSFile.FIELD_NAMES

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


def check_org_settings(org_key,setting_list) :
    settings = get_org_settings(org_key)
    return all([s in settings.keys() for s in setting_list])


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



