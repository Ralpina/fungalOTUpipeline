#!/bin/bash
#
#SBATCH --chdir=/home/yourusername/myco/
#SBATCH --job-name=cluster_denovo
#SBATCH --error=/home/yourusername/myco/errors/errorcluster.txt
#SBATCH --output=/home/yourusername/myco/errors/outputcluster.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=short
echo "hello world"

module load vsearch/2.22.1

vsearch --cluster_size ./results/notmatched_filtered.fasta --id 0.97 \
        --centroids results/denovo_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --sizeout \
        --uc results/clusters.uc


# Before assuming these de novo sequences are new OTUs, we want to make sure the quality of the original sequences is very high
# and therefore we want to inspect the chromatograms of the de novo centroids carefully
	
# Let's extract the list of de novo centroids...
grep "^>ITS" results/denovo_centroids.fasta | cut -d ";" -f 1 | sed 's/^>//' > ./results/denovosinglets.txt

# we'll use the lists to retrieve the chromatogram sequences

# now actually retrieving the sequences and copying them in the newly created folder "chrom_check"
mkdir ./results/chrom_check

singlets=($(cat ./results/denovosinglets.txt))
for singlet in "${singlets[@]}" ; do
 cp ./seqs/$singlet ./results/chrom_check
done		


# Alternatively, we can also only check the sequences of singletons:
# Let's extract the list of de novo centroids...
# grep "^>ITS" results/denovo_centroids.fasta | grep "size=1" | cut -d ";" -f 1 | sed 's/^>//' > results/denovosinglets=1.txt

# we'll use the lists to retrieve the chromatogram sequences

# now actually retrieving the sequences and copying them in the newly created folder "chrom_check_sing"
# mkdir results/chrom_check_sing

# singlets=($(cat results/denovosinglets=1.txt))
# for singlet in "${singlets[@]}" ; do
#  cp ./seqs/$singlet results/chrom_check_sing
# done		







