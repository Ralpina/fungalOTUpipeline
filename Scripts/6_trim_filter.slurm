#!/bin/bash
#
#SBATCH --chdir=/data/projects/gargiulo_popgen/myco/
#SBATCH --job-name=trim
#SBATCH --error=/data/projects/gargiulo_popgen/myco/errors/errortrim.txt
#SBATCH --output=/data/projects/gargiulo_popgen/myco/errors/outputtrim.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=short
echo "hello world"


# Now, we need to reverse-complement the reverse sequences.
# We also perform  a first trimming of the low-quality ends

# we use seqtk to reverse complement the reverse reads and to do some preliminary quality trimming on the singlets; we keep a copy of the original fastq just in case

module load seqtk/1.3
seqtk seq -r ./sing_fastq/singletsR.fastq >  ./sing_fastq/singletsReversed.fastq

# runs very quickly. Great! (sequences are reversed-complemented, quality scores are reversed)

# Seqtk uses the "modified Mott algorithm" which is the same used by phred.
# The Mott algorithm trims from both ends and has fewer issues than the sliding window used in trimmomatic. 
# averages quality scores over windows of 5 bases, and cuts when the average falls below 20 (default value)

seqtk trimfq ./sing_fastq/singletsReversed.fastq > ./sing_fastq/singletsReversedTrimmed.fastq
seqtk trimfq ./sing_fastq/singletsF.fastq > ./sing_fastq/singFTrimmed.fastq

# we can now append forward and reverse sequences to the same file

cat ./sing_fastq/singFTrimmed.fastq ./sing_fastq/singletsReversedTrimmed.fastq > ./sing_fastq/singlets.fastq



module load trimmomatic/0.39
module load vsearch/2.22.1

# for contigs
trimmomatic SE -phred33 \
        assembled_fastq/contigs.fastq \
        assembled_fastq/clean_contigs.fastq \
        LEADING:30 TRAILING:30 MINLEN:100

# for singlets
trimmomatic SE -phred33 \
        sing_fastq/singlets.fastq \
        sing_fastq/clean_singlets.fastq \
        LEADING:25 TRAILING:25 MINLEN:100

# Getting the final fastq and fasta files for Clustering and Database matching
# we keep contigs and singlets separately

# Drop the quality scores to get a fasta file for the next steps (but still keeping the fastq versions for when needed)

vsearch --fastq_filter ./assembled_fastq/clean_contigs.fastq --fastq_qmax 93 --fastaout ./assembled_fastq/clean_contigs.fasta
vsearch --fastq_filter ./sing_fastq/clean_singlets.fastq --fastq_qmax 93 --fastaout ./sing_fastq/clean_singlets.fasta
# (all sequences are kept with this last bit, it's just a file conversion)

# for the next steps
cp ./assembled_fastq/clean_contigs.fasta ./results
cp ./sing_fastq/clean_singlets.fasta ./results
cp assembled_fastq/clean_contigs.fastq ./results
cp sing_fastq/clean_singlets.fastq ./results


# let's also deal with the duplicates
# let's check how many duplicates in the singlets (works as long as there are no dots in the unique code! Otherwise change delimiter)

# we create a list of duplicate codes:
grep "^>" results/clean_singlets.fasta | cut -d "." -f 1 | cut -d "_" -f 2-10 | sort | uniq -d > duplicates.txt


# we use the tool below to generate a list of (all) singlets and their respective sequence lengths
module load seqkit/0.13.2
seqkit fx2tab --length --name --header-line  results/clean_singlets.fasta > ./results/temp.txt
cut -f 2 ./results/temp.txt | grep -v "length" > ./results/length.txt
cut -d " " -f 1 ./results/temp.txt | grep -v "#name" > ./results/name.txt
paste ./results/name.txt ./results/length.txt > ./results/sing_lengths.txt
rm ./results/temp.txt


# then we use the following loop to extract only the duplicate templates (two or more singlets for the same PCR template)
duplicates=($(cat duplicates.txt))

for duplicate in "${duplicates[@]}"; do
 grep "$duplicate" ./results/sing_lengths.txt
done > ./results/duplicate_lengths.txt

# the loop below will find the trace peak area ratio from the phred output files (this specific loop lasts for long)

singname=($(cat ./results/duplicate_lengths.txt | cut -f 1))

for name in "${singname[@]}"; do
 grep  "TRACE_PEAK_AREA_RATIO:" ./phred_out/*"$name".phd.1 | cut -d ":" -f 2
done > ./results/duplicate_peak_area.txt

# we then combine lengths and peak area ratios
paste ./results/duplicate_lengths.txt ./results/duplicate_peak_area.txt > ./results/duplicate_lengths_peaks.txt


# the script below, extracts the info about peak area ratio and length to select the best
# singlet among duplicate singlets. In particular, it takes the singlet with the smallest trace peak area ratio
# but only if the singlet is longer than 150 bp

# defining input and output files:

input_file="./results/duplicate_lengths_peaks.txt"
output_file="./results/duplicate_selected_temp.txt"

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
cut -f 1 ./results/duplicate_selected_temp.txt > ./results/selected_singlets.txt



# changing the fasta format from multi-line to single-line:

module load seqtk/1.3
seqtk seq ./results/clean_contigs.fasta > ./results/clean_contigs1.fasta
seqtk seq ./results/clean_singlets.fasta > ./results/clean_singlets1.fasta

# only renaming for contigs

cp ./results/clean_contigs1.fasta ./results/clean_contigs.fasta

# we want to remove the duplicates for the singlets
# we want to select those in the list ./results/selected_singlets.txt
# plus, we want to go back to the list "duplicates.txt", to select the singlets that had no duplicates as well

# the list of all singlets is in ./results/name.txt

duplicates=($(cat duplicates.txt))
for duplicate in "${duplicates[@]}"; do
 grep $duplicate ./results/name.txt 
done | sort  > ./results/duplicate_singlets.txt

sort ./results/name.txt

# excluding the duplicate singlets and saving only the unique ones:
sort ./results/name.txt > results/name_sorted.txt
comm -1 -3 ./results/duplicate_singlets.txt ./results/name_sorted.txt > ./results/unique_singlets.txt


# extracting the sequences of unique singlets
unique=($(cat ./results/unique_singlets.txt))
for uni in "${unique[@]}"; do
 grep -A 1 $uni ./results/clean_singlets1.fasta >> ./results/clean_unique_singlets.fasta
done


# extracting the sequences of duplicated and filtered singlets

selected=($(cat ./results/selected_singlets.txt))
for sel in "${selected[@]}"; do
 grep -A 1 $sel ./results/clean_singlets1.fasta >> ./results/clean_filt_dupl_singlets.fasta
done


# concatenating all the final singlets
cat ./results/clean_unique_singlets.fasta ./results/clean_filt_dupl_singlets.fasta > ./results/clean_filt_singlets.fasta

