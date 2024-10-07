import sys
from Bio import SeqIO
# Get lengths of all scaffolds in a fasta file (specified by varable one)
#Store results in secound variable
#Run: python get_scaffold_lengths_from_fasta.py path_to_genome scaffold_summaries/genome_name
FastaFile = open(sys.argv[1], 'r')
out_file = sys.argv[2]
with open(out_file, "a") as out_f:
    for rec in SeqIO.parse(FastaFile, 'fasta'):
        name = rec.id
        seq = rec.seq
        seqLen = len(rec)
        out_string = str(name) +","+ str(seqLen) + "\n"
        out_f.write(out_string)

FastaFile.close()