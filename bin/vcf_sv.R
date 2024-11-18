#!/usr/bin/env Rscript

# Load necessary library
library(tidyverse)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]

#genes <- args[3]

#Remove the | characters from the gene_list
#gene_list <- readLines("gene_list")
#gene_list <- gsub("\\|", "", gene_list)

# Read the sv table into a dataframe; no headers are provided
sv <- read.delim(input_file, header = F)

# Rename the columns for easier reference
sv <- sv %>% 
  rename(CHROM = V1, POS = V2, ID = V3, REF = V4, ALT = V5, END = V6, ANN = V7, LOF = V8, MEAN_COVERAGE = V9, `%VAF` = V10, SVTYPE = V11, SVLEN = V12, PHASE = V13)

# Find the maximum number of comma-separated annotations (ANN field) in the 'ANN' column
max_ann <- max(str_count(sv$ANN, ","), na.rm = T) + 1

# Split the 'ANN' column by commas into multiple columns
# Each gene/annotation will be moved into its own row
sv_ann <- sv %>%
  separate(ANN, into = paste0("ann", 1:max_ann), sep = ",") %>%
  pivot_longer(ann1:paste0("ann", max_ann), names_to = "ann_num", values_to = "ann") %>%
  filter(!is.na(ann))  # Keep only non-NA annotation rows

# Find the maximum number of '|' separated fields in the new 'ann' column (representing more details per annotation)
max_fields <- max(str_count(sv_ann$ann, "\\|"), na.rm = T) + 1

# Split the 'ann' field into multiple columns, select only relevant columns, and clean up duplicates
sv_red <- sv_ann %>%
  separate(ann, into = paste0("ann", 1:max_fields), sep = "\\|") %>%
  select(CHROM:PHASE, ann1:ann4) %>%  # Keep only first 4 annotations
  #filter(ann4 %in% gene_list) %>% # Filter for gene panel
  filter(ann3 == "HIGH" | ann3 == "MODERATE") %>%
  distinct()  # Remove duplicate rows


# Separate LOF (Loss of Function) info into individual columns
sv_red <- sv_red %>%
  separate(LOF, into = c("LOF_gene_name", "LOF_Gene_ID", "LOF_Number_of_transcripts_in_gene", "LOF_Percent_of_transcripts_affected"), sep = "\\|") %>%
  # Select relevant columns
  select(CHROM, POS, END, ann4, SVTYPE, SVLEN, `%VAF`, MEAN_COVERAGE, ann2, ann1, REF, ann3, ID, LOF_Number_of_transcripts_in_gene, LOF_Percent_of_transcripts_affected) %>%
  rename(GENES = ann4, SEQ_ONTOLOGY = ann2, EFFECT = ann3, ALT = ann1)  # Rename columns to be more descriptive

# Convert the percentage of affected transcripts to a numeric type (removing any trailing ')')
sv_red$LOF_Percent_of_transcripts_affected <- as.numeric(gsub(")", "", sv_red$LOF_Percent_of_transcripts_affected))
sv_red <- sv_red %>% 
  mutate(`%LOF_transcripts_affected` = LOF_Percent_of_transcripts_affected*100, `%VAF` = `%VAF`*100)


#Compute average MEAN_COVERAGE: 

calculate_mean <- function(numbers) {
  # Split the string into numeric values
  num_list <- as.numeric(unlist(strsplit(numbers, ",")))
  
  # Calculate the mean of the 2nd, 3rd, and 4th numbers
  mean_value <- mean(num_list[2:4])
  
  return(mean_value)
}

sv_red$MEAN_COVERAGE <- sapply(sv_red$MEAN_COVERAGE, calculate_mean)

# Make the list of genes for each event in the same row :
sv_genes <- sv_red %>%
  arrange(GENES) %>%
  group_by(ID) %>%
  summarize(
    GENES = paste(GENES, collapse = ", "),
    .groups = 'drop')

# Remove duplicated in Genes
remove_duplicates <- function(gene_col) {
  split_vector <- unlist(strsplit(gene_col, ", "))  # Split the string
  unique_values <- unique(split_vector)  # Get unique values
  result <- paste(unique_values, collapse = ", ")  # Combine back into a single string
  return(result)
}

genes_unique <- sapply(sv_genes$GENES,remove_duplicates)

sv_genes <- sv_genes %>%
  mutate(GENES = genes_unique)

sv_final <- sv_red %>%
  select(-GENES) %>%
  distinct() %>%
  left_join(sv_genes) %>%
  select(CHROM:END,GENES,SVTYPE:ID)   #Drop the LOF columns

# Replace "X" and "Y" chromosomes with numeric values for sorting (X=23, Y=24)
x <- 23
y <- 24

# Convert 'CHROM' to numeric for sorting, replacing 'X' with 23 and 'Y' with 24
sv_final$CHROM <- sv_final$CHROM %>%
  gsub("X", x, .) %>%
  gsub("Y", y, .) %>%
  gsub("chr", "", .) %>%
  as.numeric()

# Reorder the data by chromosome and position
sv_final <- sv_final %>% arrange(CHROM, POS)

# Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y chromosomes
sv_final$CHROM <- sv_final$CHROM %>%
  paste0("chr", .) %>%
  gsub(x, "X", .) %>%
  gsub(y, "Y", .)

write_tsv(sv_final, output_file)