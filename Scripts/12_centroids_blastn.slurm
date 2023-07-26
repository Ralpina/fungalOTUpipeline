module load blast/2.13.0+

split -l 20 results/denovo_centroids.fasta chunk 

for chn in chunk*; do
    echo $chn
    time blastn -remote -query $chn -db nr -num_alignments 10 -num_descriptions 10 >> results/denovo_blast.txt 
done


for chn in chunk*; do
    rm $chn
done
