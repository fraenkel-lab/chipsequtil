from StringIO import StringIO
from chipsequtil.seq import FASTAFile, FASTQFile

fasta_str = StringIO(">seq1\nACATAGGGAT\n>seq2\nTTATNTAGATA\n")
fasta_f = FASTAFile(fasta_str)
print fasta_f.headers

print "[r for r in fasta_f]", [r for r in fasta_f]
print "fasta_f['seq1']", fasta_f['seq1']
print "fasta_f.headers", fasta_f.headers
print "fasta_f.sequences", fasta_f.sequences

fastq_str = StringIO("@seq1\nACATAGGGAT\n+seq2\nY^_cccQ\JQ\n@seq2\nTTATNTAGATA\n+seq2\nY^_cJcQQJQ")
fastq_f = FASTQFile(fastq_str)
print "[r for r in fastq_f]", [r for r in fastq_f]
print "fastq_f['seq1']", fastq_f['seq1']
print "fastq_f.headers", fastq_f.headers
print "fastq_f.sequences", fastq_f.sequences
print "fastq_f.quals", fastq_f.quals
