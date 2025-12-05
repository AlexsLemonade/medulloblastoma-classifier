# J. Taroni
# 2025
#
# Extract celldex reference for use with SingleR and train SingleR model for
# use in annotating Smart-seq2 dataset
#
# USAGE: Rscript extract_ref_and_train_singler.R

# To save cache location not interactively
ExperimentHub::setExperimentHubOption("ASK", FALSE)

#### Directories and files -----------------------------------------------------

# Directory that holds processed data for single-cell experiments
single_cell_dir <- here::here("processed_data/single_cell")

# This file contains the gene identifiers that are used in the Smart-seq2
# experiment
pseudobulk_file <- fs::path(
  single_cell_dir,
  "GSE119926/GSE119926_pseudobulk_genex.tsv"
)

# Output with the model
model_file <- fs::path(
  single_cell_dir,
  "BlueprintEncode_SingleR_model.rds"
)

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
