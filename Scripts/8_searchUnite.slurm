#!/bin/bash
#
#SBATCH --chdir=/home/yourusername/myco/
#SBATCH --job-name=searchUnite
#SBATCH --error=/home/yourusername/myco/errors/errorsearchU.txt
#SBATCH --output=/home/yourusername/myco/errors/outputsearchU.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=short
echo "hello world"


module load vsearch/2.22.1

# for the contigs
vsearch --usearch_global results/clean_nochim_contigs.fasta \
        --db database/sh_general_release_dynamic_s_25.07.2023.fasta \
        --id 0.97 --strand both \
        --uc results/SH_table_contigs.uc \
        --matched results/matched_contigs.fasta \
        --notmatched results/notmatched_contigs.fasta

# for the singlets
vsearch --usearch_global results/clean_nochim_singlets.fasta \
        --db database/sh_general_release_dynamic_s_25.07.2023.fasta \
        --id 0.97 --strand both \
        --uc results/SH_table_singlets.uc \
        --matched results/matched_singlets.fasta \
        --notmatched results/notmatched_singlets.fasta
		

# We will work on the following file generated above:
# results/notmatched_contigs.fasta

# let's retrieve the list of codes for the contigs not found in UNITE at 97%
grep "^>" results/notmatched_contigs.fasta | cut -f 1 | cut -d " " -f 1 | sed 's/^>//' | cut -d "." -f 1 > ./results/notmatched_contigs.txt


# we use this list to retrieve the phd.1 files output by phred, and to re-run phd2fasta

module load phd2fasta/0.130911
mkdir results/disassembled
mkdir results/disassembled/singF
mkdir results/disassembled/singR


# retrieving files using the list 
while read string; do
files=$(find ./phred_out/ -name "*${string}*.ab1.phd.1")
for file in $files; do
 cp $file ./results/disassembled
done
done < ./results/notmatched_contigs.txt

# just for tidiness:
mv ./results/disassembled/ITS4* ./results/disassembled/singR

# running phd2fasta to concatenate all sequences in a fasta file
phd2fasta -id ./results/disassembled -os results/disassembled/singF/singF.fasta -oq results/disassembled/singF/singF.qual
phd2fasta -id ./results/disassembled/singR -os ./results/disassembled/singR/singR.fasta -oq ./results/disassembled/singR/singR.qual

# converting fasta and qual to fastq using the script "8b_fasta_to_fastq.py" (make sure it is in the myco folder)
module load python/2.7.18
python 8b_fasta_to_fastq.py

# concatenating the newly created files:
cat results/disassembled/singF/singF.fastq results/disassembled/singR/singR.fastq > results/disassembled/new_singlets.fastq

# now proceeding to quality filtering of all these new singlets:

module load trimmomatic/0.39
module load vsearch/2.22.1

trimmomatic SE -phred33 \
        ./results/disassembled/new_singlets.fastq \
        ./results/disassembled/clean_new_singlets.fastq \
        LEADING:30 TRAILING:30 MINLEN:100

# converting to fasta again

vsearch --fastq_filter ./results/disassembled/clean_new_singlets.fastq --fastq_qmax 93 --fastaout ./results/disassembled/clean_new_singlets.fasta


# Now let's deal with duplicates within these disassembled contigs.
# The script above replicates what we previously did using the script 6 on the duplicate singlets, this time on the duplicate singlets deriving from the disassembled contigs.

# Let's check how many duplicates in the singlets (the following works as long as there are no dots in the unique code! Otherwise change delimiter)

# we create a list of duplicate codes:
grep "^>" ./results/disassembled/clean_new_singlets.fasta | cut -d "." -f 1 | cut -d "_" -f 2-10 | sort | uniq -d > ./results/disassembled/duplicates.txt


# we use the tool below to generate a list of (all) singlets and their respective sequence lengths
module load seqkit/0.13.2
seqkit fx2tab --length --name --header-line  ./results/disassembled/clean_new_singlets.fasta > ./results/disassembled/temp.txt
cut -f 2 ./results/disassembled/temp.txt | grep -v "length" > ./results/disassembled/length.txt
cut -d " " -f 1 ./results/disassembled/temp.txt | grep -v "#name" > ./results/disassembled/name.txt
paste ./results/disassembled/name.txt ./results/disassembled/length.txt > ./results/disassembled/sing_lengths.txt
rm ./results/disassembled/temp.txt

# then we use the following loop to extract only the duplicate templates (two or more singlets for the same PCR template)
duplicates=($(cat ./results/disassembled/duplicates.txt))

for duplicate in "${duplicates[@]}"; do
 grep "$duplicate" ./results/disassembled/sing_lengths.txt
done > ./results/disassembled/duplicate_lengths.txt

# the loop below will find the trace peak area ratio from the phred output files (this specific loop lasts for long)

singname=($(cat ./results/disassembled/duplicate_lengths.txt | cut -f 1))

for name in "${singname[@]}"; do
 grep  "TRACE_PEAK_AREA_RATIO:" ./phred_out/*"$name".phd.1 | cut -d ":" -f 2
done > ./results/disassembled/duplicate_peak_area.txt

# we then combine lengths and peak area ratios
paste ./results/disassembled/duplicate_lengths.txt ./results/disassembled/duplicate_peak_area.txt > ./results/disassembled/duplicate_lengths_peaks.txt


# the script below extracts the info about peak area ratio and length to select the best
# singlet among duplicate singlets. In particular, it takes the singlet with the smallest trace peak area ratio
# but only if the singlet is longer than 150 bp

# defining input and output files:

input_file="./results/disassembled/duplicate_lengths_peaks.txt"
output_file="./results/disassembled/duplicate_selected_temp.txt"

# creating an empty temporary file
temp_file=$(mktemp)

# looping through each line of the input file
while IFS=$'\t' read -r name value1 value2; do
  # extracting the code from the name
  code=$(echo "$name" | cut -d "_" -f 2-20)

  # searching for rows with the same code
  same_code_rows=$(grep "$code" "$input_file")

  # if there is only one row with this code, add it to the temporary file
  if [[ $(echo "$same_code_rows" | wc -l) -eq 1 ]]; then
    echo "$same_code_rows" >> "$temp_file"
  else
    # find the row to select
    if [[ $(echo "$same_code_rows" | cut -f2 | awk '$1 > 150' | wc -l) -gt 0 ]]; then
      # select the row with the smallest value in the third column
      best_row=$(echo "$same_code_rows" | awk '$2 > 150' | sort -k3 -n | head -n1)
    else
      # select the row with the largest value in the second column
      best_row=$(echo "$same_code_rows" | awk '$2 <= 150' | sort -k2 -n | tail -n1)
      # if the values in the second columns are smaller than 150, keep the row where the value in the third column is smaller
      if [[ $(echo "$same_code_rows" | cut -f2 | awk '$1 <= 150' | wc -l) -gt 0 ]]; then
        best_row=$(echo "$same_code_rows" | awk '$2 <= 150' | sort -k3 -n | head -n1)
      fi
    fi
    # add the best row to the temporary file
    echo "$best_row" >> "$temp_file"
  fi
done < "$input_file"

# sorting the temporary file and remove duplicates
sort -u "$temp_file" > "$output_file"

# removing the temporary file
rm "$temp_file"



# let's print the list of selected singlets to a file
cut -f 1 ./results/disassembled/duplicate_selected_temp.txt > ./results/disassembled/selected_singlets.txt


# changing the fasta format from multi-line to single-line:

module load seqtk/1.3
seqtk seq ./results/disassembled/clean_new_singlets.fasta > ./results/disassembled/clean_new_singlets1.fasta

# we want to remove the duplicates for the singlets, 
# by selecting those in the list ./results/disassembled/selected_singlets.txt
# plus, we want to go back to the list "duplicates.txt", to select the singlets that had no duplicates as well

# the list of all singlets is in ./results/disassembled/name.txt

duplicates=($(cat ./results/disassembled/duplicates.txt))
for duplicate in "${duplicates[@]}"; do
 grep $duplicate ./results/disassembled/name.txt 
done | sort  > ./results/disassembled/duplicate_singlets.txt


# excluding the duplicate singlets and saving only the unique ones:
sort ./results/disassembled/name.txt > ./results/disassembled/name_sorted.txt
comm -1 -3 ./results/disassembled/duplicate_singlets.txt ./results/disassembled/name_sorted.txt > ./results/disassembled/unique_singlets.txt


# extracting the sequences of unique singlets
unique=($(cat ./results/disassembled/unique_singlets.txt))
for uni in "${unique[@]}"; do
 grep -A 1 $uni ./results/disassembled/clean_new_singlets1.fasta >> ./results/disassembled/clean_new_unique_singlets.fasta
done



# extracting the sequences of duplicated and filtered singlets

selected=($(cat ./results/disassembled/selected_singlets.txt))
for sel in "${selected[@]}"; do
 grep -A 1 $sel ./results/disassembled/clean_new_singlets1.fasta >> ./results/disassembled/clean_new_filt_dupl_singlets.fasta
done



# concatenating all the final singlets
cat ./results/disassembled/clean_new_unique_singlets.fasta ./results/disassembled/clean_new_filt_dupl_singlets.fasta > ./results/disassembled/clean_new_filt_singlets.fasta


# new file we need to search against UNITE: ./results/disassembled/clean_new_filt_singlets.fasta


# for the singlets
vsearch --usearch_global results/disassembled/clean_new_filt_singlets.fasta \
        --db database/sh_general_release_dynamic_s_25.07.2023.fasta \
        --id 0.97 --strand both \
        --uc results/disassembled/SH_table_singlets.uc \
        --matched results/disassembled/matched_singlets.fasta \
        --notmatched results/disassembled/notmatched_singlets.fasta


# now we need to concatenate the table/files generated previously for the original 

# first we concatenate the SH tables (resulting from the searches against UNITE)
# (notice that the ones not matching UNITE will be also reported and these will create duplicated accessions in the table; they will be excluded in the next steps))
cat ./results/SH_table_contigs.uc ./results/SH_table_singlets.uc ./results/disassembled/SH_table_singlets.uc > ./results/SH_table.uc

		
#let's grep SH codes and taxon names to a file:
grep -v "*" ./results/SH_table.uc | cut -f 10 | cut -d '|' -f 1,3 | sort | uniq  > ./results/SH_codes.txt


# we do not concatenate the fasta for the matching files, to avoid duplicates; they will be available as:
# results/matched_contigs.fasta 
# results/matched_singlets.fasta 
# results/disassembled/matched_singlets.fasta

# we concatenate the notmatched.fasta (we don't need the contigs because they were "disassembled" in their original singlets:
# results/notmatched_singlets.fasta
# results/disassembled/notmatched_singlets.fasta

cat ./results/notmatched_singlets.fasta ./results/disassembled/notmatched_singlets.fasta > ./results/notmatched.fasta






