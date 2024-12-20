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
result <- table_ann[
  !Annotation %in% c("intron_variant", "downstream_gene_variant", "upstream_gene_variant", "synonymous_variant") &
    grepl(paste0("(^|&)", paste(genes, collapse = "|"), "($|&)"), Gene_Name)
]
result <- result[, .(
  ID, CHROM, POS, SVTYPE, SVLEN, END, Gene_Name, REF, ALT,
  N_READS_REF, N_READS_ALT, VAF, GENOTYPE, CHR2, Annotation,
  QUAL, FILTER, Annotation_Impact, `ERRORS/WARNINGS/INFO`
)][, unique(.SD)]

# Replace "X" and "Y" chromosomes with numeric values for sorting
x <- 23
y <- 24
result[, CHROM := as.numeric(gsub("chr", "", gsub("Y", y, gsub("X", x, CHROM))))]

# Reorder by chromosome and position
setorder(result, CHROM, POS)

# Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y
result[, CHROM := gsub(y, "Y", gsub(x, "X", paste0("chr", CHROM)))]

return(result)
}


snp_arrange <- function(vcf, genes) {

  # Rename columns using data.table's setnames function
  setnames(vcf, old = grep(".*\\.GT$", colnames(vcf), value = TRUE), new = "GENOTYPE")
  setnames(vcf, old = grep(".*\\.GQ$", colnames(vcf), value = TRUE), new = "GENOTYPE_QUALITY")
  setnames(vcf, old = grep(".*\\.DP$", colnames(vcf), value = TRUE), new = "READ_DEPTH")
  setnames(vcf, old = grep(".*\\.AD$", colnames(vcf), value = TRUE), new = "ALLELE_DEPTH")
  setnames(vcf, old = grep(".*\\.AF$", colnames(vcf), value = TRUE), new = "VAF")
  
  chunk_size <- 1e6  # Number of rows per chunk
  total_rows <- nrow(vcf)

  result <- data.table()

  for (start_row in seq(1, total_rows, by = chunk_size)) {
    end_row <- min(start_row + chunk_size - 1, total_rows)
  
    # Process chunk in place
    chunk <- vcf[start_row:end_row]
  
  # Calculate maximum annotations and split
    max_ann <- max(str_count(chunk$ANN, ","), na.rm = TRUE) + 1
    ann_cols <- paste0("ann", 1:max_ann)
    chunk[, (ann_cols) := tstrsplit(ANN, ",", fixed = TRUE)]
  
    # Melt directly into the table
    chunk <- melt(chunk, measure.vars = ann_cols, value.name = "ann", na.rm = TRUE)
  
    # Split the 'ann' column
    chunk[, c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", 
              "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", 
              "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length", "CDS.pos/CDS.length", 
              "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO") := 
            tstrsplit(ann, "\\|", fixed = FALSE)]
  
    # Apply filtering and directly bind to the result
    chunk <- chunk[
      !Annotation %in% c("intron_variant", "downstream_gene_variant", "upstream_gene_variant", "synonymous_variant") &
        Gene_Name %in% genes & 
        !str_detect(Feature_ID, "XM") & 
        (Annotation_Impact %in% c("MODERATE", "HIGH")) & 
        FILTER == "PASS",
      .(ID, CHROM, POS, REF, Allele, GENOTYPE, GENOTYPE_QUALITY, READ_DEPTH, 
        ALLELE_DEPTH, VAF, QUAL, FILTER, Annotation, Annotation_Impact, Gene_Name, 
        Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVS.c, HGVS.p, 
        `cDNA.pos/cDNA.length`, `CDS.pos/CDS.length`, `AA.pos/AA.length`, Distance, 
        `ERRORS/WARNINGS/INFO`)
    ]
  
    # Append filtered rows to the result data.table
    result <- rbind(result, chunk, use.names = TRUE)
  
    gc()  # Trigger garbage collection
  }
  
  
  # Replace "X" and "Y" chromosomes with numeric values for sorting
  x <- 23
  y <- 24
  
  # Convert 'CHROM' to numeric for sorting, replacing 'X' with 23 and 'Y' with 24
  result[, CHROM := as.numeric(gsub("chr", "", gsub("X", x, gsub("Y", y, CHROM))))]
  
  # Reorder the data by chromosome and position
  setorder(result, CHROM, POS)
  
  # Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y chromosomes
  result[, CHROM := gsub(x, "X", gsub(y, "Y", paste0("chr", CHROM)))]
  
  return(result)
}

clinvar_arrange <- function(vcf, genes) {

  # Rename columns using data.table's setnames function
  setnames(vcf, old = grep(".*\\.GT$", colnames(vcf), value = TRUE), new = "GENOTYPE")
  setnames(vcf, old = grep(".*\\.GQ$", colnames(vcf), value = TRUE), new = "GENOTYPE_QUALITY")
  setnames(vcf, old = grep(".*\\.DP$", colnames(vcf), value = TRUE), new = "READ_DEPTH")
  setnames(vcf, old = grep(".*\\.AD$", colnames(vcf), value = TRUE), new = "ALLELE_DEPTH")
  setnames(vcf, old = grep(".*\\.AF$", colnames(vcf), value = TRUE), new = "VAF")
  
  chunk_size <- 20000  # Proceed in chunks to reduce memory usage
  total_rows <- nrow(vcf)
  result <- data.table()

  for (start_row in seq(1, total_rows, by = chunk_size)) {
    end_row <- min(start_row + chunk_size - 1, total_rows)
  
    # Process chunk in place
    chunk <- vcf[start_row:end_row]
  
  # Calculate maximum annotations and split
    max_ann <- max(str_count(chunk$ANN, ","), na.rm = TRUE) + 1
    ann_cols <- paste0("ann", 1:max_ann)
    chunk[, (ann_cols) := tstrsplit(ANN, ",", fixed = TRUE)]
  
    # Melt directly into the table
    chunk <- melt(chunk, measure.vars = ann_cols, value.name = "ann", na.rm = TRUE)
  
    # Split the 'ann' column
    chunk[, c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", 
              "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", 
              "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length", "CDS.pos/CDS.length", 
              "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO") := 
            tstrsplit(ann, "\\|", fixed = FALSE)]
  
    # Apply filtering and directly bind to the result
    chunk <- chunk[
      !Annotation %in% c("intron_variant", "downstream_gene_variant", "upstream_gene_variant", "synonymous_variant") &
        Gene_Name %in% genes & 
        !str_detect(Feature_ID, "XM") & 
        (Annotation_Impact %in% c("MODERATE", "HIGH")) & 
        FILTER == "PASS",
      .(ID, CHROM, POS, REF, Allele, GENOTYPE, GENOTYPE_QUALITY, READ_DEPTH, 
        ALLELE_DEPTH, VAF, QUAL, FILTER, Annotation, Annotation_Impact, Gene_Name, 
        Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVS.c, HGVS.p, 
        `cDNA.pos/cDNA.length`, `CDS.pos/CDS.length`, `AA.pos/AA.length`, Distance, 
        `ERRORS/WARNINGS/INFO`, ALLELEID, SCIDNINCL, CLNREVSTAT, ONCREVSTAT, CLNDNINCL,
        ONCDNINCL, SCIREVSTAT, CLNDN, ONCDN, SCIDN)
    ]
  
    # Append filtered rows to the result data.table
    result <- rbind(result, chunk, use.names = TRUE)
  
    gc()  # Trigger garbage collection
  }
  
  
  # Replace "X" and "Y" chromosomes with numeric values for sorting
  x <- 23
  y <- 24
  
  # Convert 'CHROM' to numeric for sorting, replacing 'X' with 23 and 'Y' with 24
  result[, CHROM := as.numeric(gsub("chr", "", gsub("X", x, gsub("Y", y, CHROM))))]
  
  # Reorder the data by chromosome and position
  setorder(result, CHROM, POS)
  
  # Convert 'CHROM' back to character format with "chr" prefix, restoring X and Y chromosomes
  result[, CHROM := gsub(x, "X", gsub(y, "Y", paste0("chr", CHROM)))]
  
  return(result)
}

remove_duplicates <- function(gene_col) {
  split_vector <- unlist(strsplit(gene_col, ", "))  # Split the string
  unique_values <- unique(split_vector)  # Get unique values
  result <- paste(unique_values, collapse = ", ")  # Combine back into a single string
  return(result)
}

spectre_arrange <- function(vcf, genes) {
  
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
    gene_pattern <- paste(genes, collapse = "|")
    result <- table_ann[str_detect(Gene_Name, gene_pattern)]
  } else {
    result <- table_ann
  }
  
  # Create list of genes in the same row, removing duplicates
  cnv_genes <- result[, .(GENES = paste(unique(Gene_Name), collapse = ", ")), by = ID]
  cnv_genes[, GENES := sapply(GENES, remove_duplicates)]
  
  # Merge gene information and select columns
  cnv_final <- merge(result[, !"Gene_Name"], cnv_genes, by = "ID", all.x = TRUE)
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
  
  # Rename 'sample_id' to 'GENOTYPE' and remove unnecessary columns

  if (any(str_detect(colnames(vcf), "reads"))) {
    setnames(vcf, old = "reads", new = "GENOTYPE")
  } else {
    sample_id <- sub("\\.wf.*", "", input_file)
    setnames(vcf, old = sample_id, new = "GENOTYPE")
  }

  vcf <- vcf[, !c("ID", "FORMAT", "INFO"), with = FALSE]
  
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

str_arrange <- function(vcf) {
  # Max values for splitting
  max_info <- max(str_count(vcf$INFO, ";"), na.rm = TRUE) + 1
  
  # Split INFO column
  info_cols <- paste0("info", seq_len(max_info))
  vcf[, (info_cols) := tstrsplit(INFO, ";", fixed = TRUE)]
  
  # Extract and rename INFO columns
  first_row <- vcf[1, ..info_cols]
  new_names <- str_extract(first_row, "^[^=]+")
  new_names[3] <- "REF_N"
  setnames(vcf, old = info_cols, new = new_names)
  
  # Process INFO columns in-place
  vcf[, (new_names) := lapply(.SD, function(col) gsub("[<>]", "", gsub("^[^=]+=", "", col))),
      .SDcols = new_names]
  
  # Rename column that contains sample ID
  if (any(str_detect(colnames(vcf), "reads"))) {
    setnames(vcf, old = "reads", new = "SAMPLE")
  } else {
    sample_id <- sub("\\.wf.*", "", input_file)
    setnames(vcf, old = sample_id, new = "SAMPLE")
  }
  
  # Split FORMAT and sample columns
  format_cols <- unique(unlist(str_split(vcf$FORMAT[1], ":", simplify = TRUE)))
  vcf[, (format_cols) := tstrsplit(SAMPLE, ":", fixed = TRUE)]
  
  # Clean ALT column
  vcf[, ALT := gsub("<STR_?([0-9]+)>", "Allele_comprised_of_\\1_repeat_units", ALT, perl = TRUE)]
  
  # Rename columns
  rename_map <- c(
    RL = "Reference_length_in_bp",
    RU = "Repeat_unit_in_the_reference_orientation",
    REPID = "REP_ID",
    VARID = "VAR_ID",
    STR_STATUS = "Repeat_expansion_status",
    STR_NORMAL_MAX = "Max_number_of_repeats_allowed_to_call_as_normal",
    STR_PATHOLOGIC_MIN = "Min_number_of_repeats_required_to_call_as_pathologic",
    SO = "Type_of_reads_that_support_the_allele",
    REPCN = "Number_of_repeat_units_spanned_by_the_allele",
    GT = "Genotype",
    REPCI = "Confidence_interval_for_REPCN",
    ADSP = "Number_of_spanning_reads_consistent_with_the_allele",
    LC = "Locus_coverage"
  )
  setnames(vcf, old = names(rename_map), new = unname(rename_map))

  
  # Remove specified columns
  cols_to_remove <- c("QUAL", "ID", "INFO", "FORMAT", "SAMPLE")
  vcf[, (cols_to_remove) := NULL]
  
  # Reorder and select columns
  setcolorder(vcf, c(
    "CHROM", "POS", "END", "REF", "ALT",
    "Repeat_unit_in_the_reference_orientation", "REP_ID",
    "Number_of_spanning_reads_consistent_with_the_allele",
    "Locus_coverage"
  ))
  
  # Process CHROM column for sorting
  vcf[, CHROM := as.integer(gsub("^chr", "", CHROM))]
  vcf[CHROM == "X", CHROM := 23]
  vcf[CHROM == "Y", CHROM := 24]
  
  # Sort by CHROM and POS
  setorder(vcf, CHROM, POS)
  
  # Restore CHROM formatting
  vcf[, CHROM := paste0("chr", CHROM)]
  vcf[CHROM == "chr23", CHROM := "chrX"]
  vcf[CHROM == "chr24", CHROM := "chrY"]
  
  return(vcf)
}

# Assuming input_file is a path to a file and vcf is a data.table
output_name <- gsub(".table", "_summary", input_file)

# Check if the data.table has any rows
if (nrow(vcf) > 0) {
  if (str_detect(input_file, "clinvar")) {
    vcf_clinvar <- clinvar_arrange(vcf, stjude_genes)  
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
  } else if (str_detect(input_file, "str")) {
    vcf_str <- str_arrange(vcf)
    fwrite(vcf_tr, paste0(output_name, ".tsv"), sep = "\t")
  } else {
    vcf_snp <- snp_arrange(vcf, stjude_genes)
    fwrite(vcf_snp, paste0(output_name, ".tsv"), sep = "\t")
  }
} else {
  message("The dataframe has no rows. No processing will be done.")
}