# for contigs:
module load vsearch/2.22.1
vsearch --uchime_ref results/clean_contigs.fasta \
        --db database/sh_general_release_dynamic_s_25.07.2023.fasta \
        --uchimealns chimeras/chimera_contigs.aln \
        --notrunclabels \
        --nonchimeras results/nochim_clean_contigs.fasta \
        --chimeras chimeras/chim_contigs.fasta \
        --borderline results/border_chim_contigs.fasta \
        --uchimeout chimeras/chim_contigs_tab.txt

# where --notrunclabels avoids truncating the header to the first space

# Borderline chimaeric sequences are sequences that have a high enough score but which are not sufficiently
# different from their closest parent.
# However, we try to be conservative and retain borderline sequences, by combining them with non-chimaeric sequences:

cat results/nochim_clean_contigs.fasta results/border_chim_contigs.fasta > results/clean_nochim_contigs.fasta


# reordering the sequences based on length first:
vsearch --sortbylength chimeras/chim_contigs.fasta \
        --output chimeras/chim_contigs_sorted.fasta \
        --fasta_width 0

# look at clustering of chimaeric sequences
vsearch --cluster_fast chimeras/chim_contigs_sorted.fasta --id 0.99 \
        --centroids chimeras/chim_contigs_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --relabel chimera_ \
        --sizeout

# for singlets:
vsearch --uchime_ref results/clean_filt_singlets.fasta \
        --db database/sh_general_release_dynamic_s_25.07.2023.fasta \
        --uchimealns chimeras/chimera_singlets.aln \
        --notrunclabels \
        --nonchimeras results/nochim_clean_singlets.fasta \
        --chimeras chimeras/chim_singlets.fasta \
        --borderline results/border_chim_singlets.fasta \
        --uchimeout chimeras/chim_singlets_tab.txt

# combine borderline with non-chimaeras
cat results/nochim_clean_singlets.fasta results/border_chim_singlets.fasta > results/clean_nochim_singlets.fasta


# reordering the sequences based on length first:
vsearch --sortbylength chimeras/chim_singlets.fasta \
        --output chimeras/chim_singlets_sorted.fasta \
        --fasta_width 0

# look at clustering of chimaeric sequences
vsearch --cluster_fast chimeras/chim_singlets_sorted.fasta --id 0.99 \
        --centroids chimeras/chim_singlets_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --relabel chimera_ \
        --sizeout

# to avoid confusion:
rm results/nochim_clean_contigs.fasta
rm results/nochim_clean_singlets.fasta
