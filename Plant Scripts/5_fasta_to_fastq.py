from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
for record in PairedFastaQualIterator(open("assembled_fastq/contigs.fasta"), open("assembled_fastq/contigs.qual")):
   print record
   
for record in PairedFastaQualIterator(open("sing_fastq/singF.fasta"), open("sing_fastq/singF.qual")):
   print record

for record in PairedFastaQualIterator(open("sing_fastq/singR.fasta"), open("sing_fastq/singR.qual")):
   print record
   


from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
handle = open("assembled_fastq/contigs.fastq", "w") #w=write
records = PairedFastaQualIterator(open("assembled_fastq/contigs.fasta"), open("assembled_fastq/contigs.qual"))
count = SeqIO.write(records, handle, "fastq")
handle.close()
print "Converted %i records" % count


handle = open("sing_fastq/singletsF.fastq", "w") #w=write
records = PairedFastaQualIterator(open("sing_fastq/singF.fasta"), open("sing_fastq/singF.qual"))
count = SeqIO.write(records, handle, "fastq")
handle.close()
print "Converted %i records" % count


handle = open("sing_fastq/singletsR.fastq", "w") #w=write
records = PairedFastaQualIterator(open("sing_fastq/singR.fasta"), open("sing_fastq/singR.qual"))
count = SeqIO.write(records, handle, "fastq")
handle.close()
print "Converted %i records" % count

