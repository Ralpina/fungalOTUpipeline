# We create directories to organise the data at the different stages of the pipeline.

mkdir plants
mkdir plants/seqs  # original sequences (scf or ab1 files), with names modified
mkdir plants/phred_out 
mkdir plants/assembled
mkdir plants/assembled_fastq
mkdir plants/chimeras
mkdir plants/results
mkdir plants/database
mkdir plants/sing
mkdir plants/singR
mkdir plants/sing_fastq
mkdir plants/scripts
mkdir plants/errors


cd plants

# IMPORTANT: At this point, we must upload the sequence files to the folder "seqs" using SFTP.
# When having many files, we can zip the original directory with the sequences to upload, and then unzip it using the command "unzip"
