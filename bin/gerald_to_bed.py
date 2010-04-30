#!/usr/bin/env python
from optparse import OptionParser
from csv import DictReader, DictWriter
import sys, os

usage = "%prog [options] <GERALD file> [<GERALD file>...]"

description = """Convert the GERALD alignment formatted files into BED format.  Input file named
<path>/<filename>.<ext> is translated into <path>/<filename>.bed unless --stdout
is specified, in which case formatted lines are written to standard output.  If
multiple input files are supplied with the --stdout option all formatted lines
are concatenated together.  Formatting only occurs for GERALD input lines that
have a valid Match Position field (i.e. successfully aligned somewhere)."""

parser = OptionParser(usage=usage, description=description)
parser.add_option('--stdout',dest='stdout',action='store_true',help='write out all formatted lines to stdout')
parser.add_option('--min-fields',dest='min_fields',action='store_true',help='only format the first three fields')
parser.add_option('--pass-only',dest='pass_only',action='store_true',help='only format lines with Y in the Pass Filtering field')

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


if __name__ == '__main__' :

	opts,args = parser.parse_args(sys.argv[1:])

	if len(args) == 0 :
		parser.print_usage()
		sys.exit(1)

	gerald_fns = args

	# step through the files
	for gerald_fn in gerald_fns :
		path,fn,fnbase,fnext = get_file_parts(gerald_fn)
		bed_lines = []


		# where to write output to
		if opts.stdout :
			f_out = sys.stdout
		else :
			f_out = open(os.path.join(path,fnbase+'.bed'),'w')

		# process input
		gerald_d = DictReader(open(gerald_fn),fieldnames=GERALDOutput.FIELD_NAMES,delimiter='\t')
		for line_d in gerald_d :
			if (opts.pass_only and line_d['filtering'] == 'Y' and line_d['match_pos'] != '') or (not opts.pass_only and line_d['match_pos'] != '') :

				outline = [line_d['match_chromo'], # chromosome
				           line_d['match_pos'], # start
				           str(int(line_d['match_pos'])+len(line_d['read'])), # end
				           line_d['read'], # read
				           '0', # score
				           '+' if line_d['match_strand'] == 'F' else '-', # strand
				           '-', # thickStart
				           '-', # thickEnd
				           '0,0,255' if line_d['match_strand'] == 'F' else '255,0,0', # itemRgb 
				          ]
				outline = '\t'.join(outline)
				f_out.write(outline+'\n')
				#bed_lines.append(bed)

		# this is the slow way
		#for line in open(gerld_fn) :
		#	grld = GERALDOutput(line)
		#	if (opts.pass_only and grld.filtering == 'Y' and grld.match_pos != '') or (not opts.pass_only and grld.match_pos != '') :
		#		bed = gerald_to_bed(grld,opts.min_fields)
		#		f_out.write(bed.output_format())
		#		#bed_lines.append(bed)

