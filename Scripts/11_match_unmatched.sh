# this script will retrieve all the sequences that were not found in UNITEv9 using the similarity threshold  >97%,
# and try to look for similar ones among the centroids previously obtained

module load vsearch/2.22.1
vsearch --usearch_global results/notmatched.fasta \
        --db results/denovo_centroids.fasta --id 0.97 \
        --fasta_width 0 \
        --uc results/matched_to_denovo.uc --strand both \
        --matched results/matched_to_denovo.fasta \
        --notmatched results/NOT_matched_to_denovo.fasta 

