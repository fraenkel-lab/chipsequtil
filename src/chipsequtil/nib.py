import math
import struct
import sys

NOMASK,MASK,HARDMASK = range(3)
def get_nib_seq(nib_fn,start,end,strand,mask=NOMASK) :
    '''Extract subsequence from .nib file like Jim Kent's nibFrag utility.

    Extract the nucleotide substring defined by the closed interval [start,end]
    from the sequence found in *nib_fn*.  *mask* parameter has the following
    possible values:

    chipsequtil.nib.NOMASK -- masked positions are not indicated (default)
    chipsequtil.nib.MASK -- masked positions are capitalized, normal bases lower case
    chipsequtil.nib.NOMASK -- masked positions are replaced with Ns
    '''

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


