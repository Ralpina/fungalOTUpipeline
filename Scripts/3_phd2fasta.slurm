#!/bin/bash
#
#SBATCH --chdir=/data/projects/gargiulo_popgen/myco/
#SBATCH --job-name=phd2fasta
#SBATCH --error=/data/projects/gargiulo_popgen/myco/errors/errorphd2fasta.txt
#SBATCH --output=/data/projects/gargiulo_popgen/myco/errors/outphd2fasta.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=all
echo "hello world"


# combining matching forward and reverse sequences in the same file using phd2fasta, i.e. preparing input files for phrap assembly	

# "phd2fasta" can produce fasta and quality files starting from phd files and can also combine forward and reverse sequences in the same file

# The unique code in the sequence names can be used to combine samples deriving from the same PCR template.
# For example, in ITS1F1_SV13S1A11.ab1.phd.1, the unique code is "SV13S1A11"

# let's extract the unique codes and use them to create a list: 

cd ./phred_out/
ls | cut -d "_" -f 2-10 | cut -d "." -f 1 | sort | uniq > ../codes.txt

# where the first cut considers "_" as a column delimiter and then cuts the second column
# the second cut considers "." as a delimiter and cuts the first column
# sort and uniq to reorder alphabetically and remove duplicates
# then creates the file codes


cd ../


# now running phd2fasta:

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

# adding sequence orientation to fasta headers:

for file in assembled/*; do
    if grep -q "^>ITS1F" $file; then
        sed -i '/^>ITS1F/ s/$/ DIRECTION: fwd/' $file
    fi
done


for file in assembled/*; do
    if grep -q "^>ITS4" $file; then
        sed -i '/^>ITS4/ s/$/ DIRECTION: rev/' $file
    fi
done
