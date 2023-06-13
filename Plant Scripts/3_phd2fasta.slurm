#!/bin/bash
#
#SBATCH --chdir=/data/projects/gargiulo_popgen/myco/plants
#SBATCH --job-name=phd2fasta
#SBATCH --error=/data/projects/gargiulo_popgen/myco/plants/errors/errorphd2fasta.txt
#SBATCH --output=/data/projects/gargiulo_popgen/myco/plants/errors/outphd2fasta.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=all
echo "hello world"


# This script will combine matching forward and reverse sequences in the same file using phd2fasta, thus preparing input files for phrap

# "phd2fasta" can produce fasta and quality files starting from phd files and can also combine forward and reverse sequences in the same file

# The unique code in the sequence names can be used to combine samples deriving from the same PCR template.
# For example, in rbcLaF_NCEA_Plant_PCR_009_H07.ab1, the unique code is "NCEA_Plant_PCR_009_H07"

# let's extract the unique codes and use them to create a list: 

cd ./phred_out/
ls | cut -d "_" -f 2-10 | cut -d "." -f 1 | sort | uniq > ../codes.txt

# where the first cut considers "_" as a column delimiter and then cuts the second column
# the second cut considers "." as a delimiter and cuts the first column
# sort and uniq to reorder alphabetically and remove duplicates
# then creates the file codes


cd ../

# combining F and R sequences in the same fasta file 

module load phd2fasta/0.130911

# Read the codes file
codes=($(cat codes.txt))

# Iterate through each code
for code in "${codes[@]}"; do
    # Find all files in the seqs_out folder that contain the code
    files=($(find phred_out -name "*$code*"))
    # Build the command
    cmd="phd2fasta ${files[@]} -os ./assembled/$code -oq ./assembled/$code.qual"
    # Execute the command
    $cmd
done

# adding orientation to fasta headers

for file in assembled/*; do
    if grep -q "^>rbcLaF" $file; then
        sed -i '/^>rbcLaF/ s/$/ DIRECTION: fwd/' $file
    fi
done


for file in assembled/*; do
    if grep -q "^>rbcLr590" $file; then
        sed -i '/^>rbcLr590/ s/$/ DIRECTION: rev/' $file
    fi
done
