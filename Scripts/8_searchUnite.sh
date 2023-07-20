module load vsearch/2.22.1

# for the contigs
vsearch --usearch_global results/clean_nochim_contigs.fasta \
        --db database/sh_general_release_dynamic_s_29.11.2022.fasta \
        --id 0.97 --strand both \
        --uc results/SH_table_contigs.uc \
        --matched results/matched_contigs.fasta \
        --notmatched results/notmatched_contigs.fasta

# for the singlets
vsearch --usearch_global results/clean_nochim_singlets.fasta \
        --db database/sh_general_release_dynamic_s_29.11.2022.fasta \
        --id 0.97 --strand both \
        --uc results/SH_table_singlets.uc \
        --matched results/matched_singlets.fasta \
        --notmatched results/notmatched_singlets.fasta
		

# we may need to change the following if using another database (or new UNITE releases)
		
#let's grep them to a file:
grep -v "*" results/SH_table_contigs.uc | cut -f 10 | cut -d '|' -f 1,3 | sort | uniq  > ./results/SH_codes_contigs.txt
grep -v "*" results/SH_table_singlets.uc | cut -f 10 | cut -d '|' -f 1,3 | sort | uniq > ./results/SH_codes_singlets.txt
# let's also concatenate them in a single file
cat ./results/SH_codes_contigs.txt ./results/SH_codes_singlets.txt | sort | uniq > ./results/SH_codes_all.txt

