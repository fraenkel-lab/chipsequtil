#!/usr/bin/env python
# nibFrag.py - a python implementation of Jim Kent's nibFrag command line utility

import struct, sys, math
from cStringIO import StringIO
from optparse import OptionParser, OptionGroup
from chipsequtil import get_file_parts

usage = '%prog [options] file.nib start end strand [outfile]'
description = """A python implementation of Jim Kent's nibFrag utility that allows outputting to
stdout.  Otherwise the functionality is identical."""
epilog="Note: not all options are currently implemented"
parser = OptionParser(usage=usage,description=description,epilog=epilog)
#parser.add_option('--output',dest='output',default=sys.stdout,help='filename to write output to [default: stdout]')
parser.add_option('--no-header',dest='no_header',action='store_true',help='only output sequence (no fasta header)')
parser.add_option('--wrap-width',dest='wrap_width',type='int',default=50,help='wrap output sequence at this number of bases, 0 indicates no wrap (sequence ends up on single line) [default: %default]')
parser.add_option('--mask-type',dest='mask_type',type='choice',choices=['NOMASK','MASK','HARDMASK'],default='NOMASK',help='how to handle masked positions, correspond to original nibFrag options --masked and --hardMasked [default: %default]')

# original nibFrag options
nibFrag_grp = OptionGroup(parser,"Original nibFrag options")
nibFrag_grp.add_option('--masked',dest='masked',action='store_true',help='use lower case characters for bases meant to be masked out')
nibFrag_grp.add_option('--hardMasked',dest='hmasked',action='store_true',help='use upper case for non masked-out and \'N\' characters for masked-out bases')
nibFrag_grp.add_option('--upper',dest='upper',action='store_true',help='use upper case characters for all bases')
nibFrag_grp.add_option('--name',dest='name',default=None,help='Use given name after \'>\' in output sequence')
parser.add_option_group(nibFrag_grp)

# nibFrag usage:
#nibFrag - Extract part of a nib file as .fa (all bases/gaps lower case by default)
#usage:
#   nibFrag [options] file.nib start end strand out.fa
#where strand is + (plus) or m (minus)
#options:
#   -masked - use lower case characters for bases meant to be masked out
#   -hardMasked - use upper case for not masked-out and 'N' characters for masked-out bases
#   -upper - use upper case characters for all bases
#   -name=name Use given name after '>' in output sequence
#   -dbHeader=db Add full database info to the header, with or without -name option
#   -tbaHeader=db Format header for compatibility with tba, takes database name as argument

NOMASK,MASK,HARDMASK = range(3)
def get_nib_seq(nib_fn,start,end,strand,mask=NOMASK) :

	nib_f = open(nib_fn,'rb')

	# first 4 bytes are a nib file signature
	#TODO - consider attempting to figure out byte order
	def_sig = 0x6BE93D3A
	sig = struct.unpack('=l',nib_f.read(4))[0]
	if def_sig != sig :
		raise Exception('Invalid nib file signature, found %s, expected %s, perhaps nib file as not created on this platform?'%(hex(sig),hex(def_sig)))

	# second 4 bytes are number of bases in sequence
	nbases = struct.unpack('=l',nib_f.read(4))[0]
	if any([nbases < c for c in [start,end]]) :
		raise Exception('Requested slice (%d,%d) not compatible with sequence of length %d in %s, aborting'%(start,end,nbases,nib_fn))

	# rest of file is sequence, with each nibble (4 bytes) being a base as follows (from http://genome.ucsc.edu/FAQ/FAQformat.html#format8) :
	#
	# 0 - T
	# 1 - C
	# 2 - A
	# 3 - G
	# 4 - N
	#
	# The most significant bit in a nibble is set if the base is masked
	trans_nuc = 'tcagn'

	# start translating the nibbles into nucleotides
	def trans_nib(nib) :
		nuc = trans_nuc[nib&7]
		mask_bit = nib & 8
		if mask in [MASK,HARDMASK] and mask_bit == 0 :
			return nuc.upper()
		if mask == HARDMASK and mask_bit == 1 :
			return 'N'
		return nuc

	# figure out how many bytes to read through
	start_byte,rem_byte = start/2,start%2
	start_out = nib_f.read(start_byte)
	seq_bytes,seq_rem_byte = int(math.ceil((end-start+rem_byte)/2.)),(end+1)%2
	seq_bytes = nib_f.read(seq_bytes+seq_rem_byte)

	# start translating the bytes
	seq = StringIO() # we use StringIO because it is more efficient than concatenating strings
	for c in seq_bytes :
		c_byte = struct.unpack('=b',c)[0]

		# higher nibble
		c_nib = (c_byte & (15<<4))>>4
		nuc = trans_nib(c_nib)
		seq.write(nuc)

		# lower nibble
		c_nib = int(c_byte) & 15
		nuc = trans_nib(c_nib)
		seq.write(nuc)

	# if the coordinate requested was not on a byte boundary, adjust
	seq_str = seq.getvalue()
	if rem_byte == 1 :
		seq_str = seq_str[1:]
	if seq_rem_byte == 1 :
		seq_str = seq_str[:-1]
	seq.close()

	return seq_str

if __name__ == '__main__' :

	opts, args = parser.parse_args(sys.argv[1:])

	if len(args) < 4 :
		parser.error('Four arguments must be supplied')

	nib_fn, strand = args[0], args[3]
	nib_path,nib_fn,nib_base,nib_ext = get_file_parts(nib_fn)
	start, stop = map(int,args[1:3])
	if stop < start :
		parser.error('Stop coordinate %d smaller than start %d'%(stop,start))

	if len(args) > 4 :
		out_f = open(args[4],'w')
	else :
		out_f = sys.stdout

	seq = get_nib_seq(nib_fn,start,stop,strand)

	out_f.write('>%s:%d-%d\n'%(nib_fn,start,stop))
	for i in xrange(0,len(seq),opts.wrap_width) :
		out_f.write(seq[i:i+opts.wrap_width]+'\n')

