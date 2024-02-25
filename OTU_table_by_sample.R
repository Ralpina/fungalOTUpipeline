setwd("C:/Users/rg53kg/OneDrive - The Royal Botanic Gardens, Kew/PacBio data trial/Preliminary - by sample/TCC samples_February 2024")

library(dplyr)
library(stringr)
library(tidyverse)

## ----------------------------------------
## LOAD UNITE9 IDENTIFICATIONS
## ----------------------------------------

library(dplyr)

# Read the data
unite <- read.delim('A11Pa.UNITE', stringsAsFactors = FALSE, header = FALSE)

unite <- unite[, c('V9', 'V4', 'V10')]
names(unite) <- c('seq', 'vsearch_perc_id', 'vsearch_identity')

# Extract numbers after "=" sign in the "seq" column (cluster sizes)
numbers <- str_extract(unite$seq, "(?<=\\=)\\d+")
# Convert to numeric, and replace non-numeric values with NA
numbers <- as.numeric(numbers)

# Add the extracted numbers as a new column to the dataframe
unite$extracted_numbers <- numbers

# Filter out rows with asterisks (*) in the "vsearch_identity" column
unite <- unite[!grepl("^\\*$", unite$vsearch_identity), ]

# Create a new dataframe to store summed abundances
summed_abundance <- unite %>%
  group_by(vsearch_identity) %>%
  summarize(abundance = sum(extracted_numbers))

# Merge the summed abundance dataframe with the original dataframe
# unite <- merge(unite, summed_abundance, by = "vsearch_identity", all.x = TRUE)

# Write the modified dataframe to a new tab-delimited file - just for tracking down things but we don't need it
# write.table(unite, file = "modified_unite_file.txt", sep = "\t", row.names = FALSE)


# Extract OTU taxonomy
unite_parts <-  strsplit(unite$vsearch_identity, '\\|')
unite$OTU <- sapply(unite_parts, function(x) paste(x[1], x[3], sep='_'))
unite$binomial <- sapply(unite_parts, '[', 1)
unite$binomial <- gsub('_', ' ', unite$binomial)


# OTU table

unite_otus <- unique(unite[,c('OTU','vsearch_identity')])
unite_taxon_parts <- strsplit(unite_otus$vsearch_identity, '\\|')
unite_taxon_parts <- strsplit(sapply(unite_taxon_parts, '[', 5), ';')
unite_taxon_parts <- matrix(unlist(unite_taxon_parts), ncol=7, byrow=TRUE)
unite_taxon_parts <- gsub('^[a-z]__', '', unite_taxon_parts)

# build taxonomy section
colnames(unite_taxon_parts) <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
unite_taxon_parts <- data.frame(unite_taxon_parts, stringsAsFactors=FALSE)

unite_otus <- cbind(unite_otus[,1,drop=FALSE], unite_taxon_parts)
unite_otus$OTU_type <- 'UNITE9'



# Add the "OTU" column to the "summed_abundance" dataframe
summed_abundance$OTU <- unite$OTU[match(summed_abundance$vsearch_identity, unite$vsearch_identity)]

# Merge "summed_abundance" with "unite_otus" by "OTU" column
unite_otus <- merge(summed_abundance, unite_otus, by = "OTU", all.x = TRUE)


write.csv(unite_otus, 'A11Pa_OTU_list.csv')



# we now transpose rows to columns and get rid of the columns we don't need

# Read the CSV file
OTU_list <- read.csv("A11Pa_OTU_list.csv", stringsAsFactors = FALSE)

# Transpose the data frame
OTU_list_transposed <- t(OTU_list)

# Remove specified rows
rows_to_remove <- c("vsearch_identity", "class", "order", "family", "genus")
OTU_list_transposed <- OTU_list_transposed[!rownames(OTU_list_transposed) %in% rows_to_remove, ]

# Use the first transposed row as the header
colnames(OTU_list_transposed) <- OTU_list_transposed[1, ]

# Remove the first row
OTU_list_transposed <- OTU_list_transposed[-1, ]

# Write the transposed data frame to a new CSV file
write.csv(OTU_list_transposed, "A11Pa_OTU_list_transposed.csv", row.names = TRUE)


# Directory containing the files
directory <- "C:/Users/rg53kg/OneDrive - The Royal Botanic Gardens, Kew/PacBio data trial/Preliminary - by sample/TCC samples_February 2024/UNITE results"

# Initialize an empty list to store OTU presence for each prefix
otu_presence <- list()

# Extract all files with the pattern "prefix_OTU_list_transposed.csv"
files <- list.files(directory, pattern = "_OTU_list_transposed.csv$", full.names = TRUE)

# Iterate over each file
for (file in files) {
    # Extract prefix from the file name
    prefix <- gsub("_OTU_list_transposed.csv", "", basename(file))
  
    # Read the second row of the CSV file to get OTU strings
    df <- read_csv(file, col_names = FALSE, skip = 1)
    otu_strings <- as.character(df[1, ])
  
    # Store OTU strings presence for the current prefix
    otu_presence[[prefix]] <- otu_strings
}

# Get unique OTU strings and sort them alphabetically
otu_strings <- sort(unique(unlist(otu_presence)))

# Initialize a list to store presence of OTU strings for each prefix
otu_matrix <- list()

# Iterate over each prefix
for (prefix in names(otu_presence)) {
    # Initialize a vector to store presence of OTU strings
    presence_vector <- numeric(length(otu_strings))
    names(presence_vector) <- otu_strings
    
    # Check if OTU string is present for the current prefix
    for (i in seq_along(otu_strings)) {
        if (otu_strings[i] %in% otu_presence[[prefix]]) {
            presence_vector[i] <- 1
        } else {
            presence_vector[i] <- 0
        }
    }
    
    # Store presence vector for the current prefix
    otu_matrix[[prefix]] <- presence_vector
}

# Create data frame from the list
df_final <- bind_cols(otu_matrix)
df_final <- as.data.frame(df_final)
rownames(df_final) <- otu_strings

# Transpose the dataframe
df_final <- t(df_final)

# Get the directory of the input files
input_dir <- dirname(files[1])

# Specify the file path for saving the CSV file in the same directory
csv_file_path <- file.path(input_dir, "otu_presence.csv")

# Print and save the dataframe
print(df_final)
write.csv(df_final, csv_file_path)

print(paste("DataFrame saved successfully to:", csv_file_path))
