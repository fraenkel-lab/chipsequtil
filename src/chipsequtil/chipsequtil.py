import os
from ConfigParser import ConfigParser


def get_file_parts(path) :
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
		return '\t'.join([str(getattr(self,d)) for d in GERALDOutput.FIELD_NAMES])+'\n'

	class FormatException(Exception) : pass


class BEDOutput :

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

	class FormatException(Exception) : pass


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
def get_org_params(org_key,addnl_configs=[]) :
  '''Returns a dict of setting/path values for a given organism that has been entered into the org_settings.py utility.'''

  
