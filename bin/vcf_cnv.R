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
if (length(args) == 3) {
    genes <- args[3]
}

# Separate the ANN annotation field into columns and keep only useful columns
cnv <- read.delim(input_file)

max_genes <- max(str_count(cnv$LOF, ","), na.rm=T) + 1

# Put each gene affected into one row and put LOF info into separated columns:
cnv <- cnv %>% 
  separate(LOF, into = paste0("gene", 1:max_genes), sep = ",") %>%
  pivot_longer(c(8:ncol(.)), names_to = "drop", values_to = "gene") %>%
  select(-drop) %>%
  filter(!is.na(gene)) %>%
  separate(gene, into = c("drop", "GENE", "N_TRANSCRIPTS", "TRANSCRIPT_LOF_(%)"), sep = "\\|") %>%
  select(-drop)

cnv$`TRANSCRIPT_LOF_(%)` <- as.numeric(gsub(")", "", cnv$`TRANSCRIPT_LOF_(%)`))

# Rename ambiguously named columns and convert LOF proportion into percentage
cnv <- cnv %>% 
mutate(`TRANSCRIPT_LOF_(%)` = `TRANSCRIPT_LOF_(%)` * 100) %>% 
rename(COPY_NUMBER_STATUS = CN)

# Reorder in function of chromosome and position
cnv$CHROM <- as.numeric(gsub("chr", "", cnv$CHROM))
cnv <- cnv %>% arrange(CHROM, POS)
cnv$CHROM <- paste0("chr", cnv$CHROM)

#Reorder the columns
cnv <- cnv %>% select(CHROM:SVTYPE,GENE,COPY_NUMBER_STATUS,N_TRANSCRIPTS,`TRANSCRIPT_LOF_(%)`)

#Keep only info on genes of interest if applicable and write table as a tsv:
if (exists("genes")) {
  gene_list <- readLines(genes)
  cnv <- cnv %>% filter(GENE %in% gene_list)
  write_tsv(cnv, output_file)
} else {
write_tsv(cnv, output_file)
}
