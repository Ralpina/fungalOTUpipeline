# Fungal OTU pipeline (for KewHPC users)

The objective of this workflow is to identify fungal Operational Taxonomic Units (OTUs) from thousands of ribosomal DNA (mainly Internal Transcribed Spacer) sequences obtained with Sanger sequencing, and to describe potentially new OTUs. Having thousands of sequences makes manual editing of chromatograms time-consuming and impractical, and therefore this workflow is recommended in that case. However, some manual editing or visual inspection of chromatograms is still recommended in the final steps of the workflow, especially when potentially describing new OTUs. 
This is by no means the best workflow; however, it is simple and can be used by people who are not familiar with scripting, provided they know the basic commands to work in a UNIX environment. 

The workflow is designed for users of the [high performance computing facility at the Royal Botanic Gardens, Kew](https://rbg-kew-bioinformatics-utils.readthedocs.io/en/latest/kewhpc/), but it could be easily adapted to other requirements, provided the main software tools are available. Most of the analyses will be run using SLURM, a job scheduler that allocates access to resources (it is good practice on KewHPC, see info [here](https://rbg-kew-bioinformatics-utils.readthedocs.io/en/latest/software/slurm/)

Software programmes and tools required (with versions used for the presented workflow):  

`phred/0.071220` (see how to obtain it [here](https://www.phrap.com/phred/))  
`phd2fasta/0.130911` (see how to obtain it [here](https://www.phrap.com/phred/))  
`phrap/1.090518` (see how to obtain it [here](https://phrap.com/))  
`python/2.7.18`  
`python/3.7.14`  
`seqtk/1.3`  
`seqkit/0.13.2`  
`vsearch/2.22.1`  
`trimmomatic/0.39`  
`fastqc/0.11.9 # (optional)`  
`openjdk # (optional)`  
`R/4.1.2 # (optional)`  
`blast/2.13.0+`  
`anaconda/2020.11`  

## Getting organised
### Sequence files naming convention
As the sequence files I inherited (almost 100K!) were not named consistently, I had to explore ways to make their analysis consistent.  
Subsequently, I built the code based on the assumption that sequence file names followed the conventions below:  
-Make sure each DNA template is indicated by a non-overlapping and unique code (i.e. a unique code MUST NOT be contained in another unique code). This unique code will become the most important identifier of your samples.  In my case, for example, "E05.ab1" is not a valid file name, as "E05" is a well name contained in other unique codes, e.g. "ITS1F_SV13S1E05.ab1" and many others;  
-The primer name must be included at the beginning of the sequence name, always followed by an underscore, as in ITS4_SV23B2A01.ab1; 
-Here we assume that two primers have been used: ITS1F for the forward sequences and ITS4 for the reverse sequences. (Whose names must be capitalised); parts of the pipeline must be changed accordingly if other primers have been used.  
-In case your unique codes are very long, avoid additional underscores within the rest of the sequence file name (there's a way to circumvent the issue of additional underscores down in the pipeline, but it's better to avoid them in the first place).  
-Avoid dots within the rest of the sequence file name (there's a way to circumvent the issue of dots down in the pipeline, but it's better to avoid them in the first place);  
-If replicates exist (i.e. sequences representing the same DNA template), add a character to distinguish them after the primer name and before the underscore and do Do not change the unique code. For example, "ITS4a_UNIQUECODE.ab1", "ITS4b_UNIQUECODE.ab1", "ITS4c_UNIQUECODE.ab1"; another option may be to add a number or letter after the unique code (as in one of the examples below), but this option shouldn't be preferred. 
-Make sure the "ab1" extension is NOT capitalised.

Following the conventions above, sequence file names will be:  
"PRIMER_UNIQUECODE*.ab1", as in:  
"ITS41_SV13S1H07.ab1"    
"ITS4_SV23B2A01_2.ab1"  
In the examples above,  
-"ITS4" is the name of the primer (the other possible primer in my sequences is ITS1F);  
-"1" after "ITS4" indicates that we can find replicate sequences from the same PCR template; their file might be named as ITS4*_ followed by a unique code ("SV13S1H07"), designated to recognise the PCR template. It is important to name replicate sequence files with the same unique code, because they will need to be analysed together;  
-likewise, the "2" in the second file indicates that we can find replicate sequences from the same PCR template; their file name will include the same unique code.   

Pro Tip: use the command `rename` to easily rename your files.  

### Running scripts and creating directories
In this workflow we will find three types of scripts:  
-slurm scripts. They can have any name and be run from any directory (provided the right working directory path is specified within the file instructions). To run a slurm script you have previously prepared following [KewHPC guidelines](https://rbg-kew-bioinformatics-utils.readthedocs.io/en/latest/software/slurm/), called "myscript.slurm":  
```sh
sbatch scripts/myscript.slurm
```
-non-slurm scripts (e.g. python or shell) have specific extension and need to made executable. For example, the scripts "myscript.py" and "myscript.sh" need to be made executable as follows:
```sh 
chmod +x myscript.py
chmod +x myscript.sh
# and then run using:
python myscript.py
./myscript.sh
```
The first script of this workflow will simply create directories to organise the data at the different stages of the pipeline. (You can always change the directory names to what works best for you).
```sh
./1_create_dirs.sh
```
The script will create the following directories:  
`myco`  
`myco/seqs`  
`myco/phred_out`  
`myco/assembled`  
`myco/assembled_fastq`  
`myco/chimeras`  
`myco/results`  
`myco/database`  
`myco/sing`  
`myco/singR`  
`myco/sing_fastq`  
`myco/scripts`  
`myco/errors`

Sequences in the format `ab1` (or `scf`) should be then uploaded to the newly created directoy `seqs`. The number of sequences uploaded can be checked using `ls seqs | wc -l`



## Basecalling using phred and first quality screening
The script "2_phred.slurm" will:  
-perform basecalling, producing output files (format "phd.1") in the folder `phred_out`. (Note that ambiguous peaks are called based on the highest peak).  
-produce a histogram with the number of high quality bases per read called "histogram_out"
-perform a first quality screening.
The quality screening will rely on the information in the phd.1 files For example:













