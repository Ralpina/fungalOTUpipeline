#!/bin/bash
#
#SBATCH --chdir=/home/yourusername/myco/
#SBATCH --job-name=blastn
#SBATCH --error=/home/yourusername/myco/errors/errorblastn.txt
#SBATCH --output=/home/yourusername/myco/errors/outputblastn.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=hmem
echo "hello world"

module load blast/2.13.0+

split -l 20 results/denovo_centroids.fasta chunk 

for chn in chunk*; do
    echo $chn
    time blastn -remote -query $chn -db nr -num_alignments 10 -num_descriptions 10 >> results/denovo_blast.txt 
done


for chn in chunk*; do
    rm $chn
done
