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


# let's concatenate the singlets and contigs files
cat results/notmatched_filtered_contigs.fasta results/notmatched_filtered_singlets.fasta > ./results/notmatched_filtered.fasta

module load vsearch/2.22.1


# clustering contigs using vsearch - only to get an idea of the relative proportion of contigs that cluster

vsearch --cluster_size results/notmatched_filtered_contigs.fasta --id 0.97 \
        --centroids results/denovo_contigs_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --sizeout

# clustering singlets using vsearch - only to get an idea of the relative proportion of singlets that cluster

vsearch --cluster_size results/notmatched_filtered_singlets.fasta --id 0.97 \
        --centroids results/denovo_singlets_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --sizeout


# clustering contigs and singlets using vsearch - gives the results we will use
	
vsearch --cluster_size ./results/notmatched_filtered.fasta --id 0.97 \
        --centroids results/denovo_centroids.fasta \
        --threads 10 \
        --strand both \
        --fasta_width 0 \
        --sizeout



# Before assuming these de novo sequences are new OTUs, we want to make sure the quality of the original sequences is very high
# and therefore we want to inspect the chromatograms of the de novo centroids carefully
	
# Let's extract the list of de novo centroids...
# ...those from singlets...	
grep "^>ITS" results/denovo_centroids.fasta | cut -d ";" -f 1 | sed 's/^>//' > denovosinglets.txt
# ...and those from contigs
grep "^>" results/denovo_centroids.fasta | grep -v "^>ITS" | cut -d ";" -f 1 | cut -d "." -f 1 | sed 's/^>//' > denovocontigs.txt	

# we'll use the lists to retrieve the chromatogram sequences

# for the contigs, we don't have info about the original sequences used to built the contig in the file name
# and we want to avoid retrieving the chromatograms for the poor quality sequences already excluded;
# so we use the phred output folder to further improve the list (the very bad sequences were filtered after phred): 		
contigs=($(cat denovocontigs.txt))		
for contig in "${contigs[@]}" ; do
 ls ./phred_out/*$contig* | cut -d "/" -f 3 | cut -d "." -f 1 >> denovocontigs2.txt 
done


# now actually retrieving the sequences and copying them in the newly created folder "chrom_check"
mkdir chrom_check

contigs2=($(cat denovocontigs2.txt))		
for contig2 in "${contigs2[@]}" ; do
 cp ./seqs/$contig2* ./chrom_check
done

singlets=($(cat denovosinglets.txt))
for singlet in "${singlets[@]}" ; do
 cp ./seqs/$singlet ./chrom_check
done		

