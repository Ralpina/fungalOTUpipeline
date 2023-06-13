#!/bin/bash
#
#SBATCH --chdir=/data/projects/gargiulo_popgen/myco/plants
#SBATCH --job-name=blastn
#SBATCH --error=/data/projects/gargiulo_popgen/myco/plants/errors/errorblastn.txt
#SBATCH --output=/data/projects/gargiulo_popgen/myco/plants/errors/outputblastn.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=hmem
echo "hello world"


# let's concatenate contigs and singlets

cat ./results/clean_contigs.fasta  ./results/clean_filt_singlets.fasta > ./results/clean_plants.fasta

module load blast/2.13.0+

split -l 20 results/clean_plants.fasta chunk 

for chn in chunk*; do
    echo $chn
    time blastn -remote -query $chn -db nr -num_alignments 10 -num_descriptions 10 >> results/denovo_blast_plants.txt 
done


for chn in chunk*; do
    rm $chn
done
