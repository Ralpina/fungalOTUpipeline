# this script will search for the centroids sequences in UNITEv9 with no similarity constraints (i.e. 0.5, which is the minimum, see vsearch manual)

module load vsearch/2.22.1
vsearch --usearch_global results/denovo_centroids.fasta \
        --db database/sh_general_release_dynamic_s_29.11.2022.fasta  \
        --id 0.5 --strand both \
        --uc results/denovo_SH_table.uc \
        --matched results/denovo_matched_SH.fasta \
        --notmatched results/denovo_NOTmatched_SH.fasta

