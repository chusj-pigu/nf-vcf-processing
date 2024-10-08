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
if(length(args) > 2) {
  genes <- args[3]
}


# Separate the ANN annotation field into columns and keep only useful columns

vcf_fusion <- read.delim(input_file)
max_ann <- max(str_count(vcf_fusion$ANN, ","), na.rm=T) + 1
vcf_fusion <- vcf_fusion %>% separate(ANN, into = paste0("ann", 1:max_ann), sep = ",") %>%
  pivot_longer(ann1:paste0("ann", max_ann), names_to = "ann_num", values_to = "ann") %>%
  filter(str_detect(ann, "gene_fusion")) %>%
  select(-ann_num) %>%
  separate(ann, into = c("Allele", "Annotation", "Annotation_impact", "Gene_name", "Gene_id", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos-length", "CDS.pos-length", "AA.pos-length", "Distance", "Message"), sep = "\\|") %>%
  separate(LOF, into = c("LOF_gene_name", "LOF_Gene_ID", "LOF_Number_of_transcripts_in_gene", "LOF_Percent_of_transcripts_affected"), sep = "\\|") %>%
  select(CHROM, POS, END, ID, Gene_name, Annotation, Annotation_impact, Allele) %>%
  rename(BREAKPOINT2=END) %>%
  distinct()

if (exists("genes")) {
  gene_fusion <- readLines(genes)
  gene_regex <- paste(gene_fusion, collapse = "|")
  vcf_fusion <- vcf_fusion %>% filter(str_detect(Gene_name, gene_regex))
  write_tsv(vcf_fusion, output_file)
} else {
  write_tsv(vcf_fusion, output_file)
}
