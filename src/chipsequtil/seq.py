from itertools import izip

def fasta_itr(f) :
    '''Returns a generator that iterates through a FASTA formatted file.
    Filename may be either a text or gzipped file.  Records are returned
    in the order they are found.'''
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
    sequences.  *fasta_fn* may be a text, gzipped file, or a file-like
    python object representing either of these.'''
    return dict(fasta_itr(f))


class FastaFile(object) :
    '''A file-like object providing information and statistics about the
    sequences in a FASTA formatted file.  Efficiently iterates through a
    text or gzipped FASTA file and provides sequential or random access to
    the records.
    
    >>> fasta_str = StringIO(">seq1\nACATAGGGAT\n>seq2\nTTATNTAGATA\n")
    >>> fasta_f = FastaFile(fasta_str)
    >>> [r for r in fasta_f]
    
    '''
    def __init__(self,f) :
        self._f = f
        self._fasta_itr = fasta_itr(f)
        self.headers = []
        self.sequences = []
        self._dict = {}

    def __get_item__(self,key) :
        return self._dict[key]

    def __set_item__(self,key,val) :
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
        return (header, seq)

