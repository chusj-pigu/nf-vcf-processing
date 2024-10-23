#!/usr/bin/env Rscript

library(tidyverse)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]

# Read the SNP table into a dataframe; no headers are provided
snp <- read.delim(input_file, header = F)

# Rename the columns for easier reference
snp <- snp %>% 
  rename(CHROM = V1, POS = V2, ID = V3, REF = V4, ALT = V5, ANN = V6, LOF = V7, GENOTYPE = V8, DEPTH_TOT = V9, `%VAF` = V10, HAPLOBLOCK = V11)

# Find the maximum number of comma-separated annotations (ANN field) in the 'ANN' column
max_ann <- max(str_count(snp$ANN, ","), na.rm = T) + 1

# Split the 'ANN' column by commas into multiple columns
# Each gene/annotation will be moved into its own row
snp_ann <- snp %>%
  select(-ID) %>%
  separate(ANN, into = paste0("ann", 1:max_ann), sep = ",") %>%
  pivot_longer(ann1:paste0("ann", max_ann), names_to = "ann_num", values_to = "ann") %>%
  filter(!is.na(ann))  # Keep only non-NA annotation rows

# Find the maximum number of '|' separated fields in the new 'ann' column (representing more details per annotation)
max_fields <- max(str_count(snp_ann$ann, "\\|"), na.rm = T) + 1

# Split the 'ann' field into multiple columns, select only relevant columns, and clean up duplicates
snp_red <- snp_ann %>%
  separate(ann, into = paste0("ann", 1:max_fields), sep = "\\|") %>%
  select(CHROM:HAPLOBLOCK, ann1:ann4) %>%  # Keep only first 4 annotations
  filter(ann3 == "HIGH" | ann3 == "MODERATE") %>% # Filter for high and medium impact
  distinct()  # Remove duplicate rows

# Separate LOF (Loss of Function) info into individual columns
snp_red <- snp_red %>%
  separate(LOF, into = c("LOF_gene_name", "LOF_Gene_ID", "LOF_Number_of_transcripts_in_gene", "%LOF_transcripts_affected"), sep = "\\|") %>%
  # Select relevant columns
  select(CHROM, POS, ann4, REF, ALT, `%VAF`, DEPTH_TOT, ann2:ann3, GENOTYPE, HAPLOBLOCK, LOF_Number_of_transcripts_in_gene, `%LOF_transcripts_affected`) %>%
  rename(GENE = ann4, SEQ_ONTOLOGY = ann2, EFFECT = ann3)  # Rename columns to be more descriptive

# Convert the percentage of affected transcripts to a numeric type (removing any trailing ')')
snp_red$`%LOF_transcripts_affected` <- as.numeric(gsub(")", "", snp_red$`%LOF_transcripts_affected`))

snp_red <- snp_red %>% 
  mutate(`%LOF_transcripts_affected` = `%LOF_transcripts_affected`*100, `%VAF` = as.numeric(`%VAF`)*100)

# Compute number of reads for alt:

snp_red <- snp_red %>%
  mutate(DEPTH_REF = ceiling(((100-`%VAF`)/100)*DEPTH_TOT)) %>%
  mutate(DEPTH_ALT = ceiling(`%VAF`*DEPTH_TOT/100)) %>%
  select(CHROM:`%VAF`,DEPTH_ALT,DEPTH_REF,DEPTH_TOT,everything())

# Replace "X" and "Y" chromosomes with numeric values for sorting (X=23, Y=24)
x <- 23
y <- 24

# Convert 'CHROM' to numeric for sorting, replacing 'X' with 23 and 'Y' with 24
snp_red$CHROM <- snp_red$CHROM %>%
  gsub("X", x, .) %>%
  gsub("Y", y, .) %>%
  gsub("chr", "", .) %>%
  as.numeric()

# Reorder the data by chromosome and position
snp_red <- snp_red %>% arrange(CHROM, POS)

# Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y chromosomes
snp_red$CHROM <- snp_red$CHROM %>%
  paste0("chr", .) %>%
  gsub(x, "X", .) %>%
  gsub(y, "Y", .)

# For now, we don't want LOF and phasing info (keeping them in case we want them in the future):

snp_red <- snp_red %>% 
  select(-`%LOF_transcripts_affected`, -LOF_Number_of_transcripts_in_gene, -HAPLOBLOCK, -GENOTYPE) %>%
  distinct()

  write_tsv(snp_red, output_file)