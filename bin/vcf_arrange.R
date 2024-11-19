#!/usr/bin/env Rscript

# Load necessary library
library(tidyverse)
library(data.table)

options(datatable.alloccol = 100)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of command line arguments are provided
 if(length(args) < 2 ){
   print("One required input missing")
   stop("Requires command line arguments: <input file> <genes for filtering CNVs and SVs> <genes for filtering SNPs>")
 }

# Assign input and output file paths
input_file <- args[1]
gene_list <- args[2]
cancer_genes <- args[3]

#Read inputs

stjude_genes <- readLines(gene_list)

sylvie_genes <- readLines(cancer_genes)

vcf <- fread(input_file)

sv_arrange <- function(vcf, genes) {

# Rename columns using patterns
colnames_to_rename <- colnames(vcf)
colnames_to_rename <- sub(".*\\.GT$", "GENOTYPE", colnames_to_rename)
colnames_to_rename <- sub(".*\\.GQ$", "GENOTYPE_QUALITY", colnames_to_rename)
colnames_to_rename <- sub(".*\\.DR$", "N_READS_REF", colnames_to_rename)
colnames_to_rename <- sub(".*\\.DV$", "N_READS_ALT", colnames_to_rename)
setnames(vcf, old = colnames(vcf), new = colnames_to_rename)

# Rename AF column to VAF
setnames(vcf, "AF", "VAF")

# Find the maximum number of comma-separated annotations (ANN field)
max_ann <- max(str_count(vcf$ANN, ","), na.rm = TRUE) + 1

# Split 'ANN' into multiple columns and reshape from wide to long format
ann_cols <- paste0("ann", 1:max_ann)
# Split the 'ANN' column into multiple columns using tstrsplit
vcf[, (ann_cols) := tstrsplit(ANN, ",", fixed = TRUE, fill = NA)]

# Reshape the data.table to long format and filter out NA values in one step
table_ann <- melt(vcf, measure.vars = ann_cols, variable.name = "ann_num", value.name = "ann")[!is.na(ann)]

# Split the 'ann' field into multiple columns and filter
ann_fields <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type",
                "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p",
                "cDNA.pos/cDNA.length", "CDS.pos/CDS.length", "AA.pos/AA.length",
                "Distance", "ERRORS/WARNINGS/INFO")
table_ann[, (ann_fields) := tstrsplit(ann, "\\|", fixed = FALSE)]

# Filter out unwanted annotations and genes, then select columns
table_red <- table_ann[
  !Annotation %in% c("intron_variant", "downstream_gene_variant", "upstream_gene_variant", "synonymous_variant") &
    Gene_Name %in% genes
]
table_red <- table_red[, .(
  ID, CHROM, POS, SVTYPE, SVLEN, END, Gene_Name, REF, ALT,
  N_READS_REF, N_READS_ALT, VAF, GENOTYPE, CHR2, Annotation,
  QUAL, FILTER, Annotation_Impact, `ERRORS/WARNINGS/INFO`
)][, unique(.SD)]

# Replace "X" and "Y" chromosomes with numeric values for sorting
x <- 23
y <- 24
table_red[, CHROM := as.numeric(gsub("chr", "", gsub("Y", y, gsub("X", x, CHROM))))]

# Reorder by chromosome and position
setorder(table_red, CHROM, POS)

# Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y
table_red[, CHROM := gsub(y, "Y", gsub(x, "X", paste0("chr", CHROM)))]

return(table_red)
}


snp_arrange <- function(vcf, genes) {

  # Rename columns using data.table's setnames function
  setnames(vcf, old = grep(".*\\.GT$", colnames(vcf), value = TRUE), new = "GENOTYPE")
  setnames(vcf, old = grep(".*\\.GQ$", colnames(vcf), value = TRUE), new = "GENOTYPE_QUALITY")
  setnames(vcf, old = grep(".*\\.DP$", colnames(vcf), value = TRUE), new = "READ_DEPTH")
  setnames(vcf, old = grep(".*\\.AD$", colnames(vcf), value = TRUE), new = "ALLELE_DEPTH")
  setnames(vcf, old = grep(".*\\.AF$", colnames(vcf), value = TRUE), new = "VAF")
  
  # Find the maximum number of comma-separated annotations (ANN field)
  max_ann <- max(str_count(vcf$ANN, ","), na.rm = TRUE) + 1
  
  # Split the 'ANN' column into multiple columns
  ann_cols <- paste0("ann", 1:max_ann)
  vcf[, (ann_cols) := tstrsplit(ANN, ",", fixed = TRUE)]

  gc()
  
  # Melt the data.table to create long format and remove NA annotations
  table_ann <- melt(vcf, measure.vars = ann_cols, value.name = "ann", na.rm = TRUE)

  gc()
  
  # Split the 'ann' column into multiple annotation fields
  table_ann[, c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", 
                "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", 
                "cDNA.pos/cDNA.length", "CDS.pos/CDS.length", "AA.pos/AA.length", 
                "Distance", "ERRORS/WARNINGS/INFO") := tstrsplit(ann, "\\|", fixed = FALSE)]

  gc()
  
  # Filter and clean the data
  table_red <- table_ann[
    !Annotation %in% c("intron_variant", "downstream_gene_variant", "upstream_gene_variant", "synonymous_variant") &
      Gene_Name %in% genes & 
      !str_detect(Feature_ID, "XM") & 
      (Annotation_Impact %in% c("MODERATE", "HIGH")),
    .(ID, CHROM, POS, REF, Allele, GENOTYPE, GENOTYPE_QUALITY, READ_DEPTH, ALLELE_DEPTH, VAF, QUAL, FILTER, 
      Annotation, Annotation_Impact, Gene_Name, Feature_Type, Feature_ID, Transcript_BioType, 
      Rank, HGVS.c, HGVS.p, `cDNA.pos/cDNA.length`, `CDS.pos/CDS.length`, `AA.pos/AA.length`, 
      Distance, `ERRORS/WARNINGS/INFO`)
  ]
  
  
  # Replace "X" and "Y" chromosomes with numeric values for sorting
  x <- 23
  y <- 24
  
  # Convert 'CHROM' to numeric for sorting, replacing 'X' with 23 and 'Y' with 24
  table_red[, CHROM := as.numeric(gsub("chr", "", gsub("X", x, gsub("Y", y, CHROM))))]
  
  # Reorder the data by chromosome and position
  setorder(table_red, CHROM, POS)
  
  # Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y chromosomes
  table_red[, CHROM := gsub(x, "X", gsub(y, "Y", paste0("chr", CHROM)))]
  
  return(table_red)
}

remove_duplicates <- function(gene_col) {
  split_vector <- unlist(strsplit(gene_col, ", "))  # Split the string
  unique_values <- unique(split_vector)  # Get unique values
  result <- paste(unique_values, collapse = ", ")  # Combine back into a single string
  return(result)
}

spectre_arrange <- function(vcf, gene_list) {
  
  # Rename columns
  setnames(vcf, old = colnames(vcf), new = gsub("^[^.]*\\.", "", colnames(vcf)))
  
  # Select relevant columns
  table <- vcf[, .(ID, CHROM, POS, END, SVLEN, SVTYPE, CN, ANN, FILTER)]
  
  # Find the maximum number of comma-separated annotations in 'ANN'
  max_ann <- max(str_count(table$ANN, ","), na.rm = TRUE) + 1
  
  # Split 'ANN' column and reshape from wide to long format
  ann_cols <- paste0("ann", 1:max_ann)
  table[, (ann_cols) := tstrsplit(ANN, ",", fixed = TRUE)]
  table_ann <- melt(table, measure.vars = ann_cols, variable.name = "ann_num", value.name = "ann")
  table_ann <- table_ann[!is.na(ann)]  # Remove rows with NA in 'ann'
  
  # Split the 'ann' field into columns and filter
  ann_fields <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", 
                  "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", 
                  "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length", "CDS.pos/CDS.length", 
                  "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO")
  table_ann[, (ann_fields) := tstrsplit(ann, "\\|", fixed = FALSE)]
  table_ann <- table_ann[
    !Annotation %in% c("intron_variant", "downstream_gene_variant", 
                       "upstream_gene_variant", "synonymous_variant") & FILTER == "PASS"
  ][, .SD, .SDcols = !c("ann_num", "Gene_ID")]
  
  # Filter based on number of rows and gene list
  if (nrow(table_ann) > 100) {
    gene_pattern <- paste(gene_list, collapse = "|")
    table_red <- table_ann[str_detect(Gene_Name, gene_pattern)]
  } else {
    table_red <- table_ann
  }
  
  # Create list of genes in the same row, removing duplicates
  cnv_genes <- table_red[, .(GENES = paste(unique(Gene_Name), collapse = ", ")), by = ID]
  cnv_genes[, GENES := sapply(GENES, remove_duplicates)]
  
  # Merge gene information and select columns
  cnv_final <- merge(table_red[, !"Gene_Name"], cnv_genes, by = "ID", all.x = TRUE)
  cnv_final <- unique(cnv_final[, .(ID, CHROM, POS, END, SVLEN, SVTYPE, CN, CHROM, FILTER, GENES, 
                                  Annotation, Annotation_Impact, `ERRORS/WARNINGS/INFO`)])
  
  # Replace "X" and "Y" chromosomes with numeric values for sorting
  x <- 23
  y <- 24
  cnv_final[, CHROM := as.numeric(gsub("chr", "", gsub("Y", y, gsub("X", x, CHROM))))]
  
  # Reorder by chromosome and position
  setorder(cnv_final, CHROM, POS)
  
  # Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y
  cnv_final[, CHROM := gsub(y, "Y", gsub(x, "X", paste0("chr", CHROM)))]
  
  return(cnv_final)
}

spectre_summary <- function(table) {
  # Select and return relevant columns
  setDT(table)
  cnv_final_red <- unique(table[, .(CHROM, POS, END, SVLEN, CN, GENES, ID)])
  return(cnv_final_red)
}


qdnaseq_arrange <- function(vcf) {
  # Convert to data.table
  setDT(vcf)
  
  # Find the maximum number of semicolon-separated values in the 'INFO' column
  max_info <- max(str_count(vcf$INFO, ";"), na.rm = TRUE) + 1
  
  # Split the 'INFO' column by semicolons into multiple columns
  info_cols <- paste0("info", 1:max_info)
  vcf[, (info_cols) := tstrsplit(INFO, ";", fixed = TRUE)]
  
  # Extract new column names from the first row of the info columns
  first_row <- vcf[1, ..info_cols]
  new_names <- sapply(first_row, function(x) str_extract(x, "^[^=]+"))
  
  # Rename the info columns with extracted names
  setnames(vcf, old = info_cols, new = new_names)
  
  # Modify each info column to contain only the text after '=' and remove angle brackets
  for (col in new_names) {
    vcf[, (col) := gsub("^[^=]+=", "", get(col))]
    vcf[, (col) := gsub("[<>]", "", get(col))]
  }
  
  # Rename 'reads' to 'GENOTYPE' and remove unnecessary columns
  setnames(vcf, "reads", "GENOTYPE")
  vcf <- vcf[, !c("ID", "FORMAT"), with = FALSE]
  
  # Replace "X" and "Y" chromosomes with numeric values for sorting (X=23, Y=24)
  x <- 23
  y <- 24
  vcf[, CHROM := gsub("X", x, CHROM)]
  vcf[, CHROM := gsub("Y", y, CHROM)]
  vcf[, CHROM := gsub("chr", "", CHROM)]
  vcf[, CHROM := as.numeric(CHROM)]
  
  # Reorder the data by chromosome and position
  setorder(vcf, CHROM, POS)
  
  # Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y chromosomes
  vcf[, CHROM := paste0("chr", CHROM)]
  vcf[, CHROM := gsub(as.character(x), "X", CHROM)]
  vcf[, CHROM := gsub(as.character(y), "Y", CHROM)]
  
  return(vcf)
}

# Assuming input_file is a path to a file and vcf is a data.table
output_name <- gsub(".table", "_summary", input_file)

# Check if the data.table has any rows
if (nrow(vcf) > 0) {
  if (str_detect(input_file, "clinvar")) {
    vcf_clinvar <- snp_arrange(vcf, stjude_genes)  
    fwrite(vcf_clinvar, paste0(output_name, ".tsv"), sep = "\t")
  } else if (str_detect(input_file, "sv")) {
    vcf_sv <- sv_arrange(vcf, stjude_genes)  
    fwrite(vcf_sv, paste0(output_name, ".tsv"), sep = "\t")
  } else if (str_detect(input_file, "spectre")) {
    vcf_spectre <- spectre_arrange(vcf,sylvie_genes)  
    vcf_spectre_red <- spectre_summary(vcf_spectre)  
    fwrite(vcf_spectre, paste0(output_name, ".tsv"), sep = "\t")
    fwrite(vcf_spectre_red, paste0(output_name, "_short.tsv"), sep = "\t")
  } else if (str_detect(input_file, "qdnaseq")) {
    vcf_qdnaseq <- qdnaseq_arrange(vcf)  
    fwrite(vcf_qdnaseq, paste0(output_name, ".tsv"), sep = "\t")
  } else {
    vcf_snp <- snp_arrange(vcf, stjude_genes)
    fwrite(vcf_snp, paste0(output_name, ".tsv"), sep = "\t")
  }
} else {
  message("The dataframe has no rows. No processing will be done.")
}