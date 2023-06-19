# let's concatenate the singlets and contigs files
cat results/notmatched_filtered_contigs.fasta results/notmatched_filtered_singlets.fasta > ./results/notmatched_filtered.fasta

module load vsearch/2.22.1


# clustering contigs using vsearch - only to get an idea of the relative proportion of contigs that cluster

vsearch --cluster_size results/notmatched_filtered_contigs.fasta --id 0.97 \
        --centroids results/denovo_contigs_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --sizeout

# clustering singlets using vsearch - only to get an idea of the relative proportion of singlets that cluster

vsearch --cluster_size results/notmatched_filtered_singlets.fasta --id 0.97 \
        --centroids results/denovo_singlets_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --sizeout


# clustering contigs and singlets using vsearch - gives the results we will use
	
vsearch --cluster_size ./results/notmatched_filtered.fasta --id 0.97 \
        --centroids results/denovo_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --sizeout
