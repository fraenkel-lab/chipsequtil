import os
from ConfigParser import ConfigParser

import chipsequtil


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
            raise GERALDOutput.FormatException('Expected %d fields in input, found %d in line: %s'%(len(GERALDOutput.FIELD_NAMES),len(line),line))

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
            raise BEDOutput.FormatException('Format requres at least 3 fields in input, found %d in line: %s'%(len(line),line))
        if len(line) > len(BEDOutput.FIELD_NAMES) :
            raise BEDOutput.FormatException('Format requres at most %d fields in input, found %d in line: %s'%(len(BEDOutput.FIELD_NAMES),len(line),line))

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


def gerald_to_bed(gerald,min_fields=False) :
    """Convert a GERALDOutput object into a BEDOutput object

    Keyword argument *min_fields* produces BED alignment with only the first three fields populated
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
                   'FDR',
                  ]

ORG_SETTINGS_FN='org_settings.cfg'
def _get_org_settings(org_key=None,addnl_configs=[]) :
    """Utility function used by get_org_settings and get_all_settings, should not be called directly"""

    config = ConfigParser()
    chipsequtil_base = os.path.split(chipsequtil.__file__)[0]
    config.read([chipsequtil_base+'/'+ORG_SETTINGS_FN,os.path.expanduser('~/.org_settings.cfg')]+addnl_configs)

    d = {}
    if org_key is None :
        for sec in config.sections() :
            d[sec] = config.items(sec)
    else :
        d = config.items(org_key)

    return dict(d)

def get_org_settings(org_key,addnl_configs=[]) :
    '''Returns a dict of setting/path values for a given organism as specified in system-wide and user's settings.'''
    return _get_org_settings(org_key,addnl_configs=[])

def get_all_settings(addnl_configs=[]) :
    '''Returns a dict of setting/path values for every organism as specified in system-wide and user's settings.'''
    return _get_org_settings(None,addnl_configs=[])


