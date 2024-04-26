#!/bin/bash
#
#SBATCH --chdir=/data/projects/gargiulo_popgen/myco/plants/assembled
#SBATCH --job-name=phrap
#SBATCH --error=/data/projects/gargiulo_popgen/myco/plants/errors/errorphrap.txt
#SBATCH --output=/data/projects/gargiulo_popgen/myco/plants/errors/outphrap.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=all
echo "hello world"


module load phrap/1.090518

# Read the codes file
codes=($(cat ../codes.txt))

# Iterate through each code
for code in "${codes[@]}"; do
    # Build the command
    cmd="phrap -bypasslevel 0 $code -ace"
    # Execute the command
    $cmd
done


echo " "


echo "Number of contigs generated from two sequences:"
cat *.ace | grep -c "AS 1 2" 
echo "Number of contigs generated from three forward or reverse sequences:"
cat *.ace | grep -c "AS 1 3" 
echo "Number of contigs generated from four forward or reverse sequences:"
cat *.ace | grep -c "AS 1 4" 
echo "Number of contigs generated from five forward or reverse sequences:"
cat *.ace | grep -c "AS 1 5"
echo "Number of contigs not generated because either forward or reverse were missing or of poor quality:"
cat *.ace | grep -c "AS 0 0" 
echo "Number of forward and reverse pairs with good quality bases but no overlap:"
cat *.ace | grep -c "AS 2" 
echo "Files with three, four, five contigs:"
cat *.ace | grep -c "AS 3" 
cat *.ace | grep -c "AS 4" 
cat *.ace | grep -c "AS 5" 
echo "done"


# notice that depending on how many times you resequenced a plate, you may have contigs generated with more than 5 forward or reverse sequences, not counted by the script above


# The script below will move contigs successfully built to the directory "assembled_fastq", whereas it will extract information about the singlets

# Loop through all ".contigs" files in the directory
for file in *.contigs; do
  # Get the prefix of the file name
  prefix=${file%.contigs}

  # Check if the corresponding ".contigs.qual" file exists
  if [ -e "${prefix}.contigs.qual" ]; then
    # Count the number of ">" characters in the file
    num_greater=$(grep -c ">" "$file")

    # If the file contains only one ">", move it to "assembled_fastq"
    if [ $num_greater -eq 1 ]; then
      mv "$file" ../assembled_fastq/
      mv "${prefix}.contigs.qual" ../assembled_fastq/
    # If the file is empty or contains two ">", add the prefix to "singlets.txt"
    elif [ $num_greater -eq 0 ] || [ $num_greater -gt 1 ]; then
      echo $prefix >> ../singlets.txt
    fi
  else
    # If the corresponding ".contigs.qual" file does not exist, add the prefix to "singlets.txt"
    echo $prefix >> ../singlets.txt
  fi
done


# we need to re-run phd2fasta for the singlets:

module load phd2fasta/0.130911

while read string; do
files=$(find ../phred_out/ -name "*${string}*.ab1.phd.1")
for file in $files; do
 cp $file ../sing
done
done < ../singlets.txt

mv ../sing/ITS3R* ../singR
phd2fasta -id ../sing -os ../sing_fastq/singF.fasta -oq ../sing_fastq/singF.qual
phd2fasta -id ../singR -os ../sing_fastq/singR.fasta -oq ../sing_fastq/singR.qual

# Understanding whether some contigs are reverse-complemented (we assume the orientation will be given by the forward primer)

# I need to extract the info from the ace file, but then need to reverse the fasta and the quality file in the folder assembled_fastq/

grep "AF ITS2F" *.ace | grep ".ab1 C" | cut -d "." -f 1 | sort | uniq > ../contigsTOreverse.txt

module load seqtk/1.3

contigs=($(cat ../contigsTOreverse.txt))
for contig in "${contigs[@]}"; do
 seqtk seq -r ../assembled_fastq/$contig.contigs > ../assembled_fastq/$contig.temp.contigs
 mv ../assembled_fastq/$contig.temp.contigs ../assembled_fastq/$contig.contigs
 seqtk seq -r ../assembled_fastq/$contig.contigs.qual > ../assembled_fastq/$contig.temp.contigs.qual
 mv ../assembled_fastq/$contig.temp.contigs.qual ../assembled_fastq/$contig.contigs.qual
done

# Let's concatenate all contigs and quality files in the same files, respectively:
cat ../assembled_fastq/*.contigs > ../assembled_fastq/contigs.fasta
cat ../assembled_fastq/*.qual > ../assembled_fastq/contigs.qual
