#!/bin/bash
#
#SBATCH --chdir=/data/projects/gargiulo_popgen/myco/
#SBATCH --job-name=hard_filt
#SBATCH --error=/data/projects/gargiulo_popgen/myco/errors/errorHardFilt.txt
#SBATCH --output=/data/projects/gargiulo_popgen/myco/errors/outputHardFilt.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=short
echo "hello world"



# this script will filter the following files 
# results/notmatched_contigs.fasta
# results/notmatched_singlets.fasta


# removing singlets with peak-area ratios > 0.15
# the command below retrieves the peak area ratios smaller than 0.15 in the original list of singlets (in awk, \n is the end of line)
awk '$3 < 0.15 {printf $1"\n"}' ./results/duplicate_lengths_peaks.txt > ./results/good_qual_singlets.txt
# the command below retrieves the list of not matched filtered singlets:
grep "^>" results/notmatched_singlets.fasta | cut -f 1 | cut -d " " -f 1 | sed 's/^>//' > ./results/notmatched_singlets.txt
# then we create a list with only the singlets that appear in both files
awk 'FNR==NR{a[$1];next}($1 in a){print}' ./results/notmatched_singlets.txt ./results/good_qual_singlets.txt > ./results/best_notmatched.txt


# we will use the ./results/best_notmatched.txt to extract the best singlets.

# going from fasta multiline to single line
module load seqtk/1.3
seqtk seq results/notmatched_contigs.fasta > results/notmatched_contigs.temp.fasta
seqtk seq results/notmatched_singlets.fasta > results/notmatched_singlets.temp.fasta

# extracting the best (not matched) singlets
while read -r line; do
 grep -A 1 $line results/notmatched_singlets.temp.fasta >> results/notmatched_singlets.temp2.fasta
done < ./results/best_notmatched.txt


# now we filter for stretches of the same nucleotides (potential indels)

grep -E 'G{9,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt
grep -E 'A{9,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt
grep -E 'T{9,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt
grep -E 'C{9,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt
grep -E 'g{9,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt
grep -E 'a{9,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt
grep -E 't{9,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt
grep -E 'c{9,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt

# filter for ambiguities (when there are more than 5 consecutive):
grep -E 'N{5,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt
grep -E 'n{5,}' -B 1 results/notmatched_contigs.temp.fasta | grep ">" | sed 's/^>//' >> results/contig_indels.txt

# -E tells grep to use extended regular expressions, which allows us to use the {} syntax to match a specific number of occurrences of a character.
# 'G{9,}' is the regular expression we're searching for. This matches any string that contains nine or more "G" characters in a row.
# with -B 1, we also grep the previous line (with the sequence name)


grep -E 'G{9,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt
grep -E 'A{9,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt
grep -E 'T{9,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt
grep -E 'C{9,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt
grep -E 'g{9,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt
grep -E 'a{9,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt
grep -E 't{9,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt
grep -E 'c{9,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt

# 48 in singlets_indels (notice that some sequences may be repeated twice because they have multiple stretches)

# filter for ambiguities (when there are more than 5 consecutive):
grep -E 'N{5,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt
grep -E 'n{5,}' -B 1 results/notmatched_singlets.temp2.fasta | grep ">" | sed 's/^>//' >> results/singlets_indels.txt


# Now we need to remove them from the contigs and singlets.

# Define the input files
indels_file="results/contig_indels.txt"
fasta_file="results/notmatched_contigs.temp.fasta"

# Define the output file
output_file="results/notmatched_filtered_contigs.fasta"

# Loop over the strings in the indels file
while read -r indel; do
  # Search for the matching string in the fasta file and remove the line and the following one
  sed -i "/$indel/{N;d;}" "$fasta_file"
done < "$indels_file"

# Copy the filtered fasta file to the output file
cp "$fasta_file" "$output_file"


# This script reads each string in contig_indels.txt and uses sed to remove the line containing the string and the following line 
# from results/notmatched_filt_contigs.temp.fasta. The modified fasta file is then saved to notmatched_filtered_contigs.fasta. 
# Note that the -i option of sed is used to modify the fasta file in-place, 
# so make sure to use a temporary file or make a backup of the original file if you need to preserve its contents.


# for the singlets 
# Define the input files
indels_file="results/singlets_indels.txt" 
fasta_file="results/notmatched_singlets.temp2.fasta"

# Define the output file
output_file="results/notmatched_filtered_singlets.fasta"

# Loop over the strings in the indels file
while read -r indel; do
  # Search for the matching string in the fasta file and remove the line and the following one
  sed -i "/$indel/{N;d;}" "$fasta_file"
done < "$indels_file"

# Copy the filtered fasta file to the output file
cp "$fasta_file" "$output_file"


# after this, we have two files for clustering:

# results/notmatched_filtered_singlets.fasta
# results/notmatched_filtered_contigs.fasta