#!/bin/bash
#
#SBATCH --chdir=/data/projects/gargiulo_popgen/myco/plants
#SBATCH --job-name=phred
#SBATCH --error=/data/projects/gargiulo_popgen/myco/plants/errors/errorphred.txt
#SBATCH --output=/data/projects/gargiulo_popgen/myco/plants/errors/out.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=all 
echo "hello world"


module load phred/0.071220

phred -id seqs -qr histogram_out -pd phred_out 

# where:
# -id indicates the directory where the files are located; 
# -qr indicates the name of the output file including a histogram with the number of high quality bases per read;
# -pd indicates the name of the directory where to save the output files in phd format (extension phd.1);

echo "phred has finished"



grep "TRIM: -1" phred_out/* | awk '{print $1}' | cut -d':' -f 1 > waste1.txt
# greps the files where the value of TRIM is minus one (by also printing their names), then pipes them to awk that only prints the first column (delimited by tab)
# then pipes them to a cut command which specifies semicolon as a delimiter and then cuts the first column with the names and path we need
# then save this list to a file called waste1.txt
grep "TRACE_PEAK_AREA_RATIO: 0.5" phred_out/* | awk '{print $1}' | cut -d':' -f 1 > waste2.txt
grep "TRACE_PEAK_AREA_RATIO: 0.4" phred_out/* | awk '{print $1}' | cut -d':' -f 1 > waste3.txt
grep "TRACE_PEAK_AREA_RATIO: 0.3" phred_out/* | awk '{print $1}' | cut -d':' -f 1 > waste4.txt
# these commands grep all the files where TRACE_PEAK_AREA_RATIO is greater than 0.3 and adds them to lists
cat waste1.txt waste2.txt waste3.txt waste4.txt > waste.txt
# the lists are then concatenated (we also keep each separate list)

# We then remove these files using the following loop:
# (be careful with rm... always make sure it is safe first, for example by first echoing the command)
# no need to add the path before $filename below, because it's already in the list waste.txt

for filename in `cat waste.txt`
  do rm $filename
done
