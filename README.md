# Fungal OTU pipeline (for KewHPC users)

The objective of this workflow is to identify fungal Operational Taxonomic Units (OTUs) from thousands of ribosomal DNA (mainly Internal Transcribed Spacer) sequences obtained with Sanger sequencing, and to describe potentially new OTUs. Having thousands of sequences makes manual editing of chromatograms time-consuming and impractical, and therefore this workflow is recommended in that case. However, some visual inspection of chromatograms is still recommended in the final steps of the workflow, especially when potentially describing new OTUs. 
This is by no means the best workflow; however, it is simple and can be used by people who are not familiar with scripting, provided they know the basic commands to work in a UNIX environment. 

The workflow is designed for users of the [high performance computing facility at the Royal Botanic Gardens, Kew](https://rbg-kew-bioinformatics-utils.readthedocs.io/en/latest/kewhpc/), but it could be easily adapted to other requirements, provided that the main software tools are available. Most of the analyses will be run using SLURM, a job scheduler that allocates access to resources (it is good practice on KewHPC, see info [here](https://rbg-kew-bioinformatics-utils.readthedocs.io/en/latest/software/slurm/)).

## Table of Contents

[Getting organised](https://github.com/Ralpina/fungalOTUpipeline#getting-organised)  

[Basecalling using phred and first quality screening](https://github.com/Ralpina/fungalOTUpipeline#basecalling-using-phred-and-first-quality-screening)  

[Preparing files for phrap assembly using phd2fasta](https://github.com/Ralpina/fungalOTUpipeline#preparing-files-for-phrap-assembly-using-phd2fasta)  

[Assembling forward and reverse sequences using phrap](https://github.com/Ralpina/fungalOTUpipeline#assembling-forward-and-reverse-files-using-phrap)  

[Merging fasta and quality files to obtain fastq files](https://github.com/Ralpina/fungalOTUpipeline#merging-fasta-and-quality-files-to-obtain-fastq-files)  

[Filtering and trimming sequences](https://github.com/Ralpina/fungalOTUpipeline#merging-fasta-and-quality-files-to-obtain-fastq-files)  

[Filtering chimaeric sequences](https://github.com/Ralpina/fungalOTUpipeline#filtering-chimaeric-sequences)

[Searching filtered sequences against the UNITE database](https://github.com/Ralpina/fungalOTUpipeline#searching-filtered-sequences-against-the-unite-database)

[Identifying sequences with no matching taxa in UNITE ("de novo" OTUs)](https://github.com/Ralpina/fungalOTUpipeline#identifying-sequences-with-no-matching-taxa-in-unite-de-novo-otus)



## Getting organised
### Sequence files naming convention

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


As the sequence files I inherited (almost 100K!) were not named consistently, I had to explore ways to make their analysis consistent. The code was then built based on the assumption that sequence file names followed the conventions below:
  - Each DNA template is indicated by a non-overlapping and unique code (i.e. a unique code MUST NOT be contained in another unique code).   
Unique code will become the most important identifier of your samples.  In my case, for example, "E05.ab1" is not a valid file name, as "E05" is a plate well name contained in other unique codes, e.g. "ITS1F_SV13S1E05.ab1" and many others;  
  - The primer name must be included at the beginning of the sequence name, always followed by an underscore, as in ITS4_SV23B2A01.ab1; 
  - Here we assume that two primers have been used: ITS1F for the forward sequences and ITS4 for the reverse sequences. (Whose names must be capitalised); parts of the pipeline must be changed accordingly if other primers have been used.  
  - In case your unique codes are very long, avoid additional underscores within the rest of the sequence file name (there's a way to circumvent the issue of additional underscores down in the pipeline, but it's better to avoid them in the first place).  
  - Avoid dots within the rest of the sequence file name (there's a way to circumvent the issue of dots down in the pipeline, but it's better to avoid them in the first place);  
  - If replicates exist (i.e. sequences representing the same DNA template), add a character to distinguish them after the primer name and before the underscore and do not change the unique code. For example, "ITS4a_UNIQUECODE.ab1", "ITS4b_UNIQUECODE.ab1", "ITS4c_UNIQUECODE.ab1". 
  - Make sure the "ab1" extension is NOT capitalised.

Following the conventions above, sequence file names will be:  
"PRIMER_UNIQUECODE*.ab1", as in:  
"ITS41_SV13S1H07.ab1"    

In the examples above,
- "ITS4" is the name of the primer (the other possible primer in my sequences is ITS1F);  
- "1" after "ITS4" indicates that we can find replicate sequences from the same PCR template; their file might be named as ITS4*_ followed by a unique code ("SV13S1H07"), designated to recognise the PCR template. It is important to name replicate sequence files with the same unique code, because they will need to be analysed together.     

Pro Tip: use the command `rename` to easily rename your files.  

### Running scripts and creating directories
In this workflow we will find three types of scripts (see "Scripts" folder):  
- slurm scripts. They can have any name and be run from any directory (provided the right working directory path is specified within the file instructions). To run a slurm script you have previously prepared following [KewHPC guidelines](https://rbg-kew-bioinformatics-utils.readthedocs.io/en/latest/software/slurm/), called "myscript.slurm":
```sh
sbatch scripts/myscript.slurm
```
- non-slurm scripts (e.g. python or shell) have specific extension and need to made executable. For example, the scripts "myscript.py" and "myscript.sh" need to be made executable as follows:
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

Sequences in the format `ab1` (or `scf`) should be then uploaded to the newly created directoy `seqs`. The number of sequences uploaded can be checked using `ls seqs | wc -l`.  
A local database can also be downloaded from the UNITE website, for example: [UNITE v.9 database](https://doi.plutof.ut.ee/doi/10.15156/BIO/2483911).



## Basecalling using phred and first quality screening
The script "2_phred.slurm" will:  
 - perform basecalling, producing output files (format "phd.1") in the folder `phred_out`. (Note that ambiguous peaks are called based on the highest peak).  
 - produce a histogram with the number of high quality bases per read called "histogram_out"  
 - perform a first quality screening.  
The quality screening will rely on the information in the phd.1 files. For example:  
```grep "TRIM: -1" phred_out/* | awk '{print $1}' | cut -d':' -f 1 > waste1.txt```  
will print the name of the files where the value of TRIM is minus one. A TRIM value equal to -1 means that the number of high quality bases is < 20, basically a file with no useful information.  
The script will also rely on the trace peak/area information in the phd.1 files, to establish whether the sequences are of sufficient quality. When peak area ratios are larger than 0.3, the sequences will be removed.  
WARNING: This quality filtering will not remove sequences with very high trace signals (and high phred scores) but potentially overlapping peaks. 

## Preparing files for phrap assembly using phd2fasta
The script "3_phd2fasta.slurm" will:  
 - extract unique codes from file names output by phred (phd.1 files);  
 - combine matching forward and reverse sequences (even if multiple ones exist due to PCR or sequencing having been repeated) in the same file, which will be named as the unique code;  
 - add information about the orientation of the sequences in the fasta headers.

## Assembling forward and reverse sequences using phrap
With the script "4_phrap_assembly.slurm":  
 - we will assemble forward and reverse sequences using phrap, recognising which ones are associated with each other by looking at unique codes;
 - phrap will generate eight files per unique code, with the following suffixes:    
```contigs```: fasta file with assembled contig(s)  
```contigs.qual```: qualities of contigs      
```5953FP.problems```: any problem encountered for this group of sequences, usually blank        
```5953FP.singlets```: fasta files with sequences not assembled  
```5953FP.contigs.qual```: qualities of contigs    
```5953FP.problems.qual```:  any problem encountered for this group of qualities, usually blank  
```5953FP.ace```: see explanation at http://bozeman.mbt.washington.edu/consed/distributions/README.29.0.txt  
```5953FP.log```: with info about the analysis.  
 - we will extract information from .ace files, about the number of contigs generated and the number of sequences from which they were generated. In particular, the first line of each ace file includes "AS <number of contigs> <total number of reads in ace file>", and can either be:  
 (1) "AS 1 2", or "AS 1 3", or "AS 1 4", depending whether there were multiple forward and reverse and one contig was generated from them    
 (2) "AS 0 0", meaning that contigs have not been produced. This happens when we have either:  
      - only one direction in the first place (only forward or only reverse sequence)  
      - bad-quality forward sequence    
      - bad-quality reverse sequence   
      - bad-quality forward and reverse  
In the cases above, the "singlets" will be saved in the file with the suffix .singlets (the reverse will not be reverse and complemented). 
 
    (3) "AS 2 2", meaning that quality for both forward and reverse was ok, but something else prevented generating a contig, probably not 
      enough overlap between the two sequences. In this case, the reverse sequences are not reverse-complemented, even if the ".contigs" file 
      is generated, and the .singlets file will be empty. However, they should be treated as singlets. It makes sense to analyse these two 
      groups of files (contigs and singlets) separately, because they will have different phred scores but also different levels of 
      "consensus".  
 - the information in the ace files will be used to move the contigs in the directory "assembled_fastq" and to create a list of singlets that will be run through phd2fasta again, to generate fasta and quality files in the folders singF and singR.  
 - for some reason, phrap reverse-complements some reads, regardless of their orientation. This information will be available in the ace files, as in:  
```AF ITS4_5952FP.ab1 C -513```  
```AF ITS1F_5952FP.ab1 U 1```  
where "C" means that the sequence has been complemented, in this case the reverse sequence. No action is needed in this case.
But in:  
```AF ITS1F_5953FP.ab1 C 1```  
```AF ITS4_5953FP.ab1 U 589```  
the script will reverse and complement the contig created (including quality file) using ```seqtk```;    
 - all contigs and singlets will be concatenated in two different fasta files (likewise quality files).  
 
 
 ## Merging fasta and quality files to obtain fastq files
The script 5_fasta_to_fastq.py simply checks that DNA and quality sequences with the same unique code are in the same order and combines them in the same fastq files for the following step. The script works in ```python/2.7``` and can be run as follows:
```sh
module load python/2.7.18
python 5_fasta_to_fastq.py
 ```
 
## Filtering and trimming sequences
The script "6_trim_filter.slurm" will:  
-reverse and complement the singlets obtained with the ITS4 (reverse primer) and trim them using ```seqtk```;  
-trim all contigs and singlets and filter out all sequences shorter than 100 nucleotides using ```trimmomatic```;  
-convert the fastq file with clean sequences to a fasta file;  
-check whether some "duplicate" templates exist, i.e. singlets obtained from the same DNA (for example too short to be assembled into a contig), and select one of them based their length (using ```seqkit``` to extract lengths) and the peak area ratio of their chromatogram. In particular, the script will keep the singlet with the smallest trace peak-area ratio, but only if the singlet is longer than 150 bp;  
 -produce the following files:  
 ```./results/clean_filt_singlets.fasta```  
 ```./results/clean_contigs.fasta```  
WARNING: This quality filtering will not remove sequences with very high trace signals (and high phred scores) but potentially overlapping peaks.

## Filtering chimaeric sequences
The script "7_chimaera_filter.sh" will search for chimaeric sequences in contigs and singlets against the [UNITE v.9 database](https://doi.plutof.ut.ee/doi/10.15156/BIO/2483911) using vsearch. At the end of the script, potentially chimaeric sequences will be filtered out and two files will be produced:  
 ```results/clean_nochim_singlets.fasta```  
  ```results/clean_nochim_contigs.fasta```


## Searching filtered sequences against the UNITE database
The script "8_searchUnite.sh" will:
- search all filtered contigs and singlets against [UNITE v.9](https://doi.plutof.ut.ee/doi/10.15156/BIO/2483911) (using the algorithm usearch_global in vsearch), using a similarity threshold of 97% (notice that this threshold can be changed based on your needs);  
- print lists of taxa with corresponding SH codes, for contigs, singlets and combined contigs/singlets (with extention .txt);  
- the script will also produce the following files:  
   ```results/SH_table_contigs.uc```: full table of contigs in the format explained in the [vsearch manual](https://vcru.wisc.edu/simonlab/bioinformatics/programs/vsearch/vsearch_manual.pdf);    
   ```results/SH_table_singlets.uc```: full table of singlets, as above;  
   ```results/matched_contigs.fasta```: sequences of the contigs for which a match in UNITE was found (over 97% similarity);   
   ```results/matched_singlets.fasta```: sequences of the singlets for which a match in UNITE was found (over 97% similarity);    
   ```results/notmatched_contigs.fasta```: sequences of the contigs for which a match in UNITE was not found (they may match with a lower similarity threshold);      
   ```results/notmatched_singlets.fasta```: sequences of the singlets for which a match in UNITE was not found (they may match with a lower similarity threshold).      

At this point, we can have an idea of the taxonomic composition of our dataset (which will depend on the taxonomic composition of UNITEv9!).
 
 
## Identifying sequences with no matching taxa in UNITE ("de novo" OTUs)
We can now work on the sequences not found in UNITE when using 97% as a similarity threshold. These sequences will be stored in:  
```results/notmatched_contigs.fasta```  
```results/notmatched_singlets.fasta```  

### Hard-filtering sequences for clustering
The script "9_hard_filt.slurm" will:    
- extract singlets with peak-area ratios < 0.15: the rationale for this filter is to get rid of all sequences for which basecalling was uncertain (those with trace peak-area ratios > 0.15), even if in small portions of the sequences. Trace peak-area ratios are retrieved from the file ./results/duplicate_lengths_peaks.txt, created in one of the previous steps;  
- convert fasta files from multi-line to single-line, search for and remove stretches of more than nine identical nucleotides (potential indels) which may have disrupted basecalling (from both contigs and singlets);    
- filter out sequences with more than five consecutive Ns from both contigs and singlets. The script will NOT look for additional ambiguities, as phred and phrap are not expected to have produced any, with the options used above;  
- produce two filtered files to be used in the subsequent clustering:  
 ```results/notmatched_filtered_singlets.fasta```   
 ```results/notmatched_filtered_contigs.fasta```  

WARNING: This quality filtering will not remove sequences with very high trace signals (and high phred scores) but potentially overlapping peaks. A MANUAL STEP WITH VISUAL INSPECTION OF CHROMATOGRAMS IS RECOMMENDED AT THIS STAGE OR AFTER CLUSTERING. 

 ### Clustering sequences to obtain centroids
The script "10_denovo_centroids.sh" will cluster sequences in vsearch based on abundance (cluster_size), using an identity threshold = 97% (note that you can adjust this percentage based on your needs). The script will output a fasta file, ```denovo_centroids.fasta```, including all the centroid sequences and the relative size of each cluster. 
WARNING: If you have not inspected chromatogram sequences of the de novo centroids, A VISUAL INSPECTION OF CHROMATOGRAMS IS RECOMMENDED NOW.

 
### Trying to assign all sequences to clusters 
The script "11_match_unmatched.sh" will try to assign all the sequences that were not found in UNITEv9 using a similarity threshold = 97%
to the centroids previously obtained. In particular, it will:  
- concatenate the files ```./results/notmatched_contigs.fasta``` and ```./results/notmatched_singlets.fasta``` in the file ```./results/notmatched.fasta``` (notice that these files derive from the script 8 and are not filtered);  
- use ```vsearch``` to compare sequences to the centroids in the file ```./results/denovo_centroids.fasta```   
- produce three files, including the list of the sequences matching to centroids (```./results/matched_to_denovo.uc```), and fasta files with sequences matching (```./results/matched_to_denovo.fasta```) or not matching to centroids (```./results/NOT_matched_to_denovo.fasta```) using a 97% similarity threshold.
 
 
### Identifying "de novo" centroid sequences
Here, two potential and alternative identification methods will be described: (1) blasting the centroid sequences against the NCBI database; and (2) searching the centroid sequences against UNITE using no similarity constraints. 
 
#### Method 1: Blasting centroid sequences against the global database  
The script "12_centroids_blastn.slurm" will:  
- carry out a remote nucleotide blast using the blast binary ```blast/2.13.0+```. Notice that this is performed by splitting the input sequence file in subsets, to decrease the computational effort;  
- produce an output called ```./results/denovo_blast_summary.txt``` with the 10 most likely blast results for each query sequence.
The output can then be used by the script "13_centroids_blast_summary.py", that will build a table with the most likely blast result for each query, its taxonomic information, and the various blast metrics.
 
#### Method 2: Searching centroid sequences against UNITE with no similarity constraints
The script "14_searchUniteFree.sh" will search the centroids sequences against UNITEv9 with no similarity constraints (i.e. 0.5, which is the minimum, see the [vsearch manual](https://vcru.wisc.edu/simonlab/bioinformatics/programs/vsearch/vsearch_manual.pdf)).
The script will produce three output files:   
- ```./results/denovo_SH_table.uc```: with the search results in the format explained in the [vsearch manual](https://vcru.wisc.edu/simonlab/bioinformatics/programs/vsearch/vsearch_manual.pdf), and including the de novo centroid name with the cluster size in column 9, the matching sequence in UNITE in column 10, and the percentage of similarity in column 4;  
- ```./results/denovo_matched_SH.fasta```: including the de novo centroid sequences with a match in UNITE;  
- ```./results/denovo_NOTmatched_SH.fasta```: including the de novo centroid sequences without a match in UNITE.


### Assigning ecological guilds to "de novo" sequences
The script 15_funguild.sh will use [FUNGuild](https://github.com/UMNFuN/FUNGuild) to assign ecological guilds to the putative de novo sequences.
...TO BE CONTINUED...




 
### Acknowledgements
Scripts number 5, 12 and 13 are associated with [van der Linde et al. 2018](https://www.nature.com/articles/s41586-018-0189-9).
Roberta Gargiulo's work is funded by Defra.

### References
- Nguyen NH, Song Z, Bates ST, Branco S, et al. 2016. FUNGuild: an open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology 20, 241-248.
- van der Linde S, Suz LM, Orme CDL, et al. 2018 Environment and host as large-scale controls of ectomycorrhizal fungi. Nature 558, 243–248 (2018). https://doi.org/10.1038/s41586-018-0189-9








