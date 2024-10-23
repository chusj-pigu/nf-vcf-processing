#!/usr/bin/env Rscript

# Load necessary library
library(tidyverse)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of command line arguments are provided
 if(length(args) < 2){
   print("Please include input and output file paths as arguments!")
   stop("Requires command line arguments: <input_file> <output_file>")
 }

# Assign input and output file paths
input_file <- args[1]
output_file <- args[2]

# Separate the ANN annotation field into columns and keep only useful columns
cnv <- read.delim(input_file, header = F)

cnv <- cnv %>% 
  rename(CHROM = V1, POS = V2, ID = V3, REF = V4, ALT = V5, END = V6, SVLEN = V7, SVTYPE = V8, CN = V9, LOF = V10)

max_genes <- max(str_count(cnv$LOF, ","), na.rm=T) + 1

# Put each gene affected into one row and put LOF info into separated columns:
cnv_ann <- cnv %>% 
  separate(LOF, into = paste0("gene", 1:max_genes), sep = ",") %>%
  pivot_longer(c(10:ncol(.)), names_to = "drop", values_to = "gene") %>%
  select(-drop) %>%
  filter(!is.na(gene)) %>%
  separate(gene, into = c("drop", "GENE", "N_TRANSCRIPTS", "TRANSCRIPT_LOF_(%)"), sep = "\\|") %>%
  select(-drop)

cnv_ann$`TRANSCRIPT_LOF_(%)` <- as.numeric(gsub(")", "", cnv_ann$`TRANSCRIPT_LOF_(%)`))

# Rename ambiguously named columns and convert LOF proportion into percentage
cnv_red <- cnv_ann %>% 
  mutate(`TRANSCRIPT_LOF_(%)` = `TRANSCRIPT_LOF_(%)` * 100) %>% 
  rename(COPY_NUMBER_STATUS = CN) %>%
  select(CHROM,POS,END:GENE,ID)     # For now, we don't want LOF info

# Make the list of genes in the same row :
cnv_genes <- cnv_red %>%
  arrange(GENE) %>%
  group_by(ID) %>%
  summarize(
    GENES = paste(GENE, collapse = ", "),
    .groups = 'drop')

cnv_final <- cnv_red %>%
  select(-GENE) %>%
  distinct() %>%
  left_join(cnv_genes) %>%
  select(CHROM:COPY_NUMBER_STATUS,GENES,ID)

# Reorder in function of chromosome and position
cnv_final$CHROM <- as.numeric(gsub("chr", "", cnv_final$CHROM))
cnv_final <- cnv_final %>% arrange(CHROM, POS)
cnv_final$CHROM <- paste0("chr", cnv_final$CHROM)

#Keep only info on genes of interest if applicable and write table as a tsv:
write_tsv(cnv_final, output_file)