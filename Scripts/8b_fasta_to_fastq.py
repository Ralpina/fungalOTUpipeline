from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
for record in PairedFastaQualIterator(open("results/disassembled/singF/singF.fasta"), open("results/disassembled/singF/singF.qual")):
   print record
   
for record in PairedFastaQualIterator(open("results/disassembled/singR/singR.fasta"), open("results/disassembled/singR/singR.qual")):
   print record

from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
handle = open("results/disassembled/singF/singF.fastq", "w") #w=write
records = PairedFastaQualIterator(open("results/disassembled/singF/singF.fasta"), open("results/disassembled/singF/singF.qual"))
count = SeqIO.write(records, handle, "fastq")
handle.close()
print "Converted %i records" % count


handle = open("results/disassembled/singR/singR.fastq", "w") #w=write
records = PairedFastaQualIterator(open("results/disassembled/singR/singR.fasta"), open("results/disassembled/singR/singR.qual"))
count = SeqIO.write(records, handle, "fastq")
handle.close()
print "Converted %i records" % count

