# Assuming we have FUNGuild.py in the directory FUNGuild, we need the following files:

# results of the search of sequences against UNITE with 0.97 similarity threshold:
cp ./results/SH_table.uc ./FUNGuild/

# results of the search of de novo sequences against UNITE with 0.5 similarity threshold:
cp results/denovo_SH_table.uc ./FUNGuild/


cd FUNGuild/
grep -v "^N" SH_table.uc | cut -f 9-10 > SH.uc
grep -v "^N" denovo_SH_table.uc | cut -f 9-10 > denovo.uc 
grep -v "^N" denovo_SH_table.uc | cut -f 4 > similar.uc 

# we add the header (make sure separators are the same as for the following lines, usually tabs)
# most importantly, we need to add "taxonomy" as the header of the taxonomic information (case-sensitive!), in the last column
sed -i $'1 i x\ttaxonomy' SH.uc
sed -i $'1 i x\ttaxonomy' denovo.uc
# header for the similarity file
sed -i '1i similarity' similar.uc

module load python/3.7.14
module load anaconda/2020.11 ## used instead of pandas on KewHPC
# preparing the taxonomic table

python FUNGuild.py taxa -otu SH.uc -format tsv -classifier unite
python FUNGuild.py taxa -otu denovo.uc -format tsv -classifier unite
# these will produce the two files SH.taxa.uc and denovo.taxa.uc
# assigning guilds
python FUNGuild.py guild -taxa SH.taxa.uc                                      
python FUNGuild.py guild -taxa denovo.taxa.uc 
# these will produce the two files: SH.taxa.guilds.txt and denovo.taxa.guilds.txt


# now let's prepare tables that also include SH codes
cut -f 2 SH.uc | cut -d "|" -f 3 > SH.txt
cut -f 2 denovo.uc | cut -d "|" -f 3 > SH_denovo.txt
paste SH.taxa.guilds.txt SH.txt > guilds.SH.txt
paste denovo.taxa.guilds.txt SH_denovo.txt similar.uc > denovo.guilds.SH.txt

# we want to extract all the accessions which are assigned with some probability to a given taxon
# we do it for the three different files, to be able to tell apart de novo from those found in UNITE
grep "Probable" guilds.SH.txt | grep "Ectomycorrhizal" > Ectoguilds.txt
grep "Probable" denovo.guilds.SH.txt | grep "Ectomycorrhizal" > Ectodenovoguilds.txt


# to count how many SHs:
echo "Total number of de novo Ectomycorrhizal sequences (centroids)"
wc -l Ectodenovoguilds.txt 
echo "Number of different SHs in all the de novo Ectomycorrhizal sequences (centroids)"
cut -f 18 Ectodenovoguilds.txt | sort | uniq | wc -l
echo "Number of different SHs in all Ectomycorrhizal sequences with 97% similarity to UNITE SHs"
cut -f 18 Ectoguilds.txt | sort | uniq | wc -l

# to extract the non-ecto sequences, we do the opposite as we did above:
grep "Probable" guilds.SH.txt | grep -v "Ectomycorrhizal" > non_ectoguilds.txt
grep "Probable" denovo.guilds.SH.txt | grep -v "Ectomycorrhizal" > non_ectodenovoguilds.txt

echo "Number of different SHs in all the de novo non-ectomycorrhizal sequences (centroids)"
cut -f 18 non_ectodenovoguilds.txt | sort | uniq | wc -l
echo "Number of different SHs in all non-ectomycorrhizal sequences with 97% similarity to UNITE SHs"
cut -f 18 non_ectoguilds.txt | sort | uniq | wc -l
