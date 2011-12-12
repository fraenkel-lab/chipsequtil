from chipsequtil import get_org_settings, BEDFile
from chipsequtil.nib import NibDB
from pprint import pprint

genome_dir = get_org_settings('mm9')['genome_dir']
db = NibDB(nib_dirs=[genome_dir])
fasta_headers,seqs = db.get_fasta_from_bed('shuffled_peaks.bed')

pprint(seqs[:10])
