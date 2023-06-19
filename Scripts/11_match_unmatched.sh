# this script will retrieve all the sequences that were not found in UNITEv9 using a similarity threshold  = 97%,
# and try to look for similar ones among the centroids previously obtained

# let's first concatenate the singlets and contigs fasta file, to avoid getting too many output files:
cat results/notmatched_contigs.fasta results/notmatched_singlets.fasta > ./results/notmatched.fasta


module load vsearch/2.22.1
vsearch --usearch_global results/notmatched.fasta \
        --db results/denovo_centroids.fasta --id 0.97 \
        --fasta_width 0 \
        --uc results/matched_to_denovo.uc --strand both \
        --matched results/matched_to_denovo.fasta \
        --notmatched results/NOT_matched_to_denovo.fasta 

