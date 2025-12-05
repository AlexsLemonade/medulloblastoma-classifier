# J. Taroni
# 2025
#
# Extract celldex reference for use with SingleR and train SingleR model and
# extract ENSG identifiers for oligodendrocyte markers from PanglaoDB
# for use in annotating Smart-seq2 dataset
#
#
# USAGE: Rscript set_up_cell_annotation_references.R

#### Directories and files -----------------------------------------------------

# Directory that holds processed data for single-cell experiments
single_cell_dir <- here::here("processed_data/single_cell")

# This file contains the gene identifiers that are used in the Smart-seq2
# experiment
pseudobulk_file <- fs::path(
  single_cell_dir,
  "GSE119926/GSE119926_pseudobulk_genex.tsv"
)

# PanglaoDB file
panglaodb_file <- here::here("data/PanglaoDB/PanglaoDB_markers_2020-03-27.tsv")

# Output with the SingleR model
model_file <- fs::path(
  single_cell_dir,
  "BlueprintEncode_SingleR_model.rds"
)

# Output of GSEABase::GeneSet with PanglaoDB markers for oligodendrocytes
oligodendrocyte_out_file <- fs::path(
  single_cell_dir,
  "PanglaoDB_oligodendrocytes_ENSG.rds"
)

#### Annotation setup ----------------------------------------------------------

# Source gene identifier conversion function
source(here::here("utils/convert_gene_names.R"))

# To save cache location not interactively
ExperimentHub::setExperimentHubOption("ASK", FALSE)

#### Extract reference and train SingleR model ---------------------------------

# Use the BlueprintEncodeData as a reference
ref <- celldex::BlueprintEncodeData(ensembl = TRUE)

# Grab the gene identifiers from the pseudobulk data
gene_ids <- readr::read_tsv(pseudobulk_file) |>
  dplyr::pull(gene)

# Train a SingleR model on BlueprintEncode data
singler_model <- SingleR::trainSingleR(
  ref = ref,
  labels = ref$label.fine,
  genes = "de",
  restrict = gene_ids
)

# Write model to file to use for classification downstream
readr::write_rds(
  singler_model,
  file = model_file,
  compress = "gz"
)

#### PanglaoDB -----------------------------------------------------------------

# Read in and subset to human markers
panglaodb_df <- readr::read_tsv(panglaodb_file) |>
  dplyr::filter(stringr::str_detect(species, "Hs"))

# Subset to oligodendrocytes markers
oligodendrocyte_markers <- panglaodb_df |>
  dplyr::filter(`cell type` == "Oligodendrocytes")

# Convert from gene symbols to Ensembl gene identifiers
oligodendrocyte_markers <- convert_gene_names(
  genex_df <- oligodendrocyte_markers,
  gene_column_before = "official gene symbol",
  gene_column_after = "gene",
  map_from = "SYMBOL",
  map_to = "ENSEMBL"
) |>
  dplyr::select(cell_type = `cell type`,
                gene)

# In the format that's needed for AUCell
oligodendrocyte_gene_set <- GSEABase::GeneSet(
  oligodendrocyte_markers$gene,
  setName = "panglaodb_oligodendrocyte",
  geneIdType = GSEABase::ENSEMBLIdentifier()
)

# Write to file
readr::write_rds(
  oligodendrocyte_gene_set,
  file = oligodendrocyte_out_file,
  compress = "gz"
)
