from itertools import izip

def fasta_itr(f) :
    '''Returns a generator that iterates through a FASTA formatted file.
    *f* may be either a text or gzipped file, or a file-like python object
    representing either of these.  Records are returned in the order they
    are found.'''
    if isinstance(f,str) :
        f = open(f)

    # check for magic number 1f 8b indicating gzip file, I dunno, just cuz
    if f.read(2) == "\x1f\x8b" : f = gzip.GzipFile(f)
    else : f.seek(0)

    curr_header, curr_seq = None, None
    for r in f :
        if r.startswith('>') :
            if curr_header is not None :
                yield (curr_header, curr_seq)
            curr_header = r[1:].strip()
            curr_seq = ''
        else :
            curr_seq += r.strip()
    # return the last record
    yield (curr_header,curr_seq)

def fasta_to_dict(f) :
    '''Returns a dictionary whose keys are FASTA headers and values are
    sequences.  *f* may be a text, gzipped file, or a file-like
    python object representing either of these.'''
    return dict(fasta_itr(f))


def fastq_itr(f) :
    '''Returns a generator that iterates through a FASTQ formatted file.
    *f* may be either a text or gzipped file, or a file-like python object
    representing either of these.  Records are returned in the order they
    are found.'''
    if isinstance(f,str) :
        f = open(f)

    # check for magic number 1f 8b indicating gzip file, I dunno, just cuz
    if f.read(2) == "\x1f\x8b" : f = gzip.GzipFile(f)
    else : f.seek(0)

    SEQ, QUAL = 0,1
    in_region = SEQ
    curr_header, curr_seq, curr_qual = None, None, None
    for r in f :
        if r.startswith('@') :
            if curr_header is not None :
                yield (curr_header, (curr_seq, curr_qual))
            curr_header = r[1:].strip()
            curr_seq = ''
            curr_qual = ''
            in_region = SEQ
        elif r.startswith('+') :
            in_region = QUAL
        else :
            curr_field = r.strip()
            if in_region == SEQ :
                curr_seq += curr_field
            elif in_region == QUAL :
                curr_qual += curr_field

    # return the last record
    yield (curr_header,(curr_seq,curr_qual))

def fastq_to_dict(f) :
    '''Returns a dictionary whose keys are FASTQ headers and values are
    sequences.  *f* may be a text, gzipped file, or a file-like
    python object representing either of these.'''
    return dict(fastq_itr(f))


class FASTAFile(object) :
    '''A file-like object providing information and statistics about the
    sequences in a FASTA formatted file.  Efficiently iterates through a
    text or gzipped FASTA file and provides sequential or random access to
    the records.  Instances store header and sequence data as they are read
    
      >>> fasta_str = StringIO(">seq1\\nACATAGGGAT\\n>seq2\\nTTATNTAGATA\\n")
      >>> fasta_f = FASTAFile(fasta_str)
      >>> [r for r in fasta_f]
      [('seq1', 'ACATAGGGAT'), ('seq2', 'TTATNTAGATA')]
      >>> fasta_f['seq1']
      ACATAGGGAT
      >>> fasta_f.headers
      ['seq1', 'seq2']
      >>> fasta_f.sequences
      ['ACATAGGGAT', 'TTATNTAGATA']

    Instances have the following members:

    **headers**
      list of FASTA headers in original order

    **sequences**
      list of FASTA sequences in original order

    .. NOTE::
       The members **headers** and **sequences** are not available until the
       the FASTA records have been iterated once

    '''

    def __init__(self,f) :
        self._f = f
        self._fasta_itr = fasta_itr(f)
        self.headers = []
        self.sequences = []
        self._dict = {}

    def __getitem__(self,key) :
        return self._dict[key]

    def __setitem__(self,key,val) :
        self._dict[key] = val

    def next(self) :
        '''Returns next FASTA record in the file as (header, sequence) tuple.'''

        if self._fasta_itr is None :
            self._fasta_itr = izip(self.headers,self.sequences)

        try :
            header, seq = self._fasta_itr.next()
        except StopIteration, e :
            self._fasta_itr = None
            self._f = None
            raise e

        if self._f is not None : 
            # this means we're not done reading through the file yet
            self.headers.append(header)
            self.sequences.append(seq)
            self._dict[header] = seq

        return header, seq

    def __iter__(self) :
        return self


class FASTQFile(object) :
    '''A file-like object providing information and statistics about the
    sequences in a FASTQ formatted file.  Efficiently iterates through a
    text or gzipped FASTQ file and provides sequential or random access to
    the records.  Instances store header and sequence data as they are read
    
      >>> fastq_str = StringIO("@seq1\\nACATAGGGAT\\n+seq2\\nY^_cccQ\JQ\\n
      @seq2\\nTTATNTAGATA\\n+seq2\\nY^_cJcQQJQ")
      >>> fastq_f = FASTQFile(fastq_str)
      >>> [r for r in fastq_f]
      [('seq1', 'ACATAGGGAT'), ('seq2', 'TTATNTAGATA')]
      >>> fastq_f['seq1']
      ACATAGGGAT
      >>> fastq_f.headers
      ['seq1', 'seq2']
      >>> fastq_f.sequences
      ['ACATAGGGAT', 'TTATNTAGATA']

    Instances have the following members:

    **headers**
      list of FASTQ headers in original order

    **sequences**
      list of FASTQ sequences in original order

    **quals**
      list of FASTQ quality scores in original order

    .. NOTE::
       The members **headers** and **sequences** are not available until the
       the FASTA records have been iterated once

    '''

    def __init__(self,f) :
        self._f = f
        self._fastq_itr = fastq_itr(f)
        self.headers = []
        self.sequences = []
        self.quals = []
        self._dict = {}

    def __getitem__(self,key) :
        return self._dict[key]

    def __setitem__(self,key,val) :
        self._dict[key] = val

    def next(self) :
        '''Returns next FASTA record in the file as (header, sequence) tuple.'''

        if self._fastq_itr is None :
            self._fastq_itr = izip(self.headers,self.sequences)

        try :
            header, (seq, qual) = self._fastq_itr.next()
        except StopIteration, e :
            self._fastq_itr = None
            self._f = None
            raise e

        if self._f is not None : 
            # this means we're not done reading through the file yet
            self.headers.append(header)
            self.sequences.append(seq)
            self.quals.append(qual)
            self._dict[header] = (seq, qual)

        return header, seq

    def __iter__(self) :
        return self

