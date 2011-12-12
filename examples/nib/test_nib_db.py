from chipsequtil import get_org_settings, BEDFile
from chipsequtil.nib import NibDB
from pprint import pprint

# see `org_settings.py -h` for more info on get_org_settings(<organism>) function
genome_dir = get_org_settings('mm9')['genome_dir']

# NibDB is an interface to a collection of nib files, typically corresponding
# to chromosomes of a genome

# example with only one nib file
print 'NibDB with a single nib file'
db = NibDB(nib_fns=[genome_dir+'/chr1.nib'])

print 'NibDB info:'
pprint(dict(db.db_info))

# get a fasta record for some sequence
print 'Example fasta record: chr1:1e8-1e8+100'
print db.get_fasta('chr1',1e8,1e8+100)

# get just the sequence
print 'Same example, only sequence:'
print db.get_seq('chr1',1e8,1e8+100)
print


# example with a directory of nib files
print 'NibDB with a directory of nib files'
db = NibDB(nib_dirs=[genome_dir])

# get a fasta record for some sequence
print 'Example fasta record: chr1:1e8-1e8+100'
print db.get_fasta('chr1',1e8,1e8+100)

print 'Example fasta record: chr1:1e8-1e8+100'
print db.get_fasta('chr2',1e8,1e8+100)

print 'Example fasta record: chr1:1e8-1e8+100'
print db.get_fasta('chrX',1e8,1e8+100)


# example of fetching all sequences from a bed file
fasta_headers,seqs = db.get_fasta_from_bed('shuffled_peaks.bed')

print 'Num. peaks:',len(open('shuffled_peaks.bed').readlines())
pprint(seqs[:10])
