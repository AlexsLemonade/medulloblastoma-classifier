# Gather metadata from each data source
#
# Steven Foltz
# November 2022

data_dir <- here::here("data")
processed_data_dir <- here::here("processed_data")
processed_pseudobulk_data_dir <- here::here(processed_data_dir, "pseudobulk")

dir.create(processed_data_dir, showWarnings = FALSE)
dir.create(processed_pseudobulk_data_dir, showWarnings = FALSE)

# input file names
GSE124814_metadata_input_filename <- here::here(data_dir, "GSE124814", "GSE124814_sample_descriptions.xlsx")
GSE164677_metadata_input_filename <- here::here(data_dir, "GSE164677", "GSE164677_series_matrix.txt.gz")
openpbta_metadata_input_filename <- here::here(data_dir, "OpenPBTA", "pbta-histologies.tsv")
sj_metadata_input_filename <- here::here(data_dir, "stjudecloud", "SAMPLE_INFO.txt")
GSE119926_metadata_input_filename <- here::here(data_dir, "GSE119926", "GSE119926_series_matrix.txt.gz")
GSE155446_metadata_input_filename <- here::here(data_dir, "GSE155446", "GSE155446_human_cell_metadata.csv.gz")

# output file names
bulk_metadata_output_filename <- here::here(processed_data_dir,
                                            "bulk_metadata.tsv")
pseudobulk_metadata_output_filename <- here::here(processed_pseudobulk_data_dir,
                                                  "pseudobulk_metadata.tsv")

################################################################################
# functions
################################################################################

clean_mb_subgroups <- function(df){

  df <- df |>
    dplyr::mutate(subgroup = dplyr::case_when(subgroup %in% c("E", "Group 3", "Group3", "Group3_alpha", "Group3_beta", "Group3_gamma", "MB_GRP3", "GP3") ~ "G3",
                                              subgroup %in% c("C", "D", "Group 4", "Group4", "Group4_alpha", "Group4_beta", "Group4_gamma", "MB_GRP4", "GP4") ~ "G4",
                                              subgroup %in% c("NORM", "n/a (NORM)") ~ "Normal",
                                              subgroup %in% c("B", "MB_SHH", "SHH_alpha", "SHH_beta", "SHH_delta", "SHH_gamma", "SHH-infant", "SHH-adult") ~ "SHH",
                                              subgroup %in% c("A", "MB_WNT", "WNT_alpha", "WNT_beta") ~ "WNT",
                                              subgroup %in% c("GP3/4") ~ NA,
                                              .default = subgroup))

  return(df)

}

################################################################################
# GSE124814
################################################################################

# Weishaupt et al. harmonized metadata for GSE124814:
#   Mapped original subgroups (subgroup_supplied_original) to consensus subgroups (subgroup_supplied_renamed)
#   Identified duplicate samples (merged duplicate samples together by averaging expression values)
#   Limited to just MB and normal cerebellum tissues
#   Citation: Weishaupt, H. et al. Batch-normalization of cerebellar and medulloblastoma gene expression datasets utilizing empirically defined negative control genes. Bioinformatics 35, 3357–3364 (2019)

# Read in GSE124814 metadata

GSE124814_samples_df_column_names <- c("sample_name",
                                       "title",
                                       "CEL_file",
                                       "source_name",
                                       "organism",
                                       "molecule",
                                       "label",
                                       "chip_name",
                                       "age",
                                       "sex", # was gender
                                       "histology",
                                       "brain_region",
                                       "death",
                                       "metastatic_stage",
                                       "follow_up",
                                       "subgroup_supplied_original",
                                       "subgroup_supplied_renamed",
                                       "subgroup_relabeled",
                                       "description",
                                       "description2")

# when reading in the xlsx file,
# 1. use the column names above
# 2. skip final column
# 3. ignore first two rows
GSE124814_metadata <- readxl::read_xlsx(GSE124814_metadata_input_filename,
                                        col_names = GSE124814_samples_df_column_names,
                                        col_types = "text",
                                        skip = 2) |>
  tidyr::separate(title, # separate title into experiment accession and sample id
                  into = c("experiment_accession", "id"),
                  sep = "-",
                  extra = "merge") |> # split on the first hyphen and keep rest intact
  # filter out E-MTAB-292 per data accessibility
  dplyr::filter(experiment_accession != "EMTAB292") |>
  # GSE124814 merged duplicate samples together by averaging expression values
  # and this is indicated in the description column with the word "average"
  dplyr::mutate(is_duplicate = stringr::str_detect(description, "average"),
                # "NA" subgroup corresponds to "Normal"
                subgroup_supplied_renamed = ifelse(subgroup_supplied_renamed == "NA",
                                                   "Normal",
                                                   subgroup_supplied_renamed),
                # treat "Unknown" subgroup as NA
                subgroup_supplied_renamed = dplyr::na_if(x = subgroup_supplied_renamed,
                                                         y = "Unknown"),
                # isolate the sample accession as its own column ("reanalysis of SAMPLE_ACCESSION (EXPERIMENT_ACCESSION)")
                sample_accession = stringr::word(description, 3) ) |>
  dplyr::select(sample_accession,
                subgroup = subgroup_supplied_renamed,
                study = experiment_accession,
                is_duplicate) |>
  dplyr::mutate(platform = "Array") |>
  clean_mb_subgroups()

################################################################################
# GSE164677
################################################################################

# Get the GEO series metadata file and transform it
GSE164677_metadata <- GEOquery::getGEO(filename = GSE164677_metadata_input_filename) |>
  as.data.frame() |>
  dplyr::select(sample_accession = geo_accession,
                subgroup = medulloblastoma.subgroup.ch1) |>
  dplyr::mutate(study = "GSE164677",
                is_duplicate = FALSE,
                platform = "RNA-seq") |>
  clean_mb_subgroups()

################################################################################
# OpenPBTA (MB)
################################################################################

openpbta_mb_metadata <- readr::read_tsv(file = openpbta_metadata_input_filename,
                                        col_types = "c") |>
  dplyr::filter(experimental_strategy == "RNA-Seq",
                short_histology == "Medulloblastoma") |>
  tidyr::separate(molecular_subtype, # separate molecular_subtype into molecular and subgroup
                  into = c("molecular", "subgroup"),
                  sep = ", ",
                  extra = "merge") |>
  dplyr::mutate(subgroup = dplyr::na_if(x = subgroup,
                                        y = "To be classified")) |>
  dplyr::arrange(Kids_First_Participant_ID, Kids_First_Biospecimen_ID) |> # patient ID, sample ID
  dplyr::mutate(is_duplicate = duplicated(Kids_First_Participant_ID)) |> # marks 2+ instance of patient ID
  dplyr::rename("sample_accession" = "Kids_First_Biospecimen_ID") |>
  dplyr::mutate(study = "OpenPBTA",
                platform = "RNA-seq") |>
  dplyr::select(sample_accession,
                subgroup,
                study,
                is_duplicate,
                platform) |>
  clean_mb_subgroups()

################################################################################
# OpenPBTA (LGG)
################################################################################

openpbta_lgg_metadata <- readr::read_tsv(file = openpbta_metadata_input_filename,
                                         col_types = "c") |>
  dplyr::filter(experimental_strategy == "RNA-Seq",
                pathology_diagnosis == "Low-grade glioma/astrocytoma (WHO grade I/II)",
                short_histology == "LGAT") |>
  dplyr::mutate(subgroup = "LGG") |>
  dplyr::arrange(Kids_First_Participant_ID, Kids_First_Biospecimen_ID) |> # patient ID, sample ID
  dplyr::mutate(is_duplicate = duplicated(Kids_First_Participant_ID)) |> # marks 2+ instance of patient ID
  dplyr::rename("sample_accession" = "Kids_First_Biospecimen_ID") |>
  dplyr::mutate(study = "OpenPBTA",
                platform = "RNA-seq") |>
  dplyr::select(sample_accession,
                subgroup,
                study,
                is_duplicate,
                platform) # do not clean_mb_subgroups()

################################################################################
# St. Jude
################################################################################

sj_metadata <- readr::read_tsv(file = sj_metadata_input_filename,
                               col_types = "c") |>
  dplyr::filter(lubridate::mdy(sj_embargo_date) |>
                  lubridate::year() < 2023) |> # keep samples with embargo ending before 2023
  dplyr::arrange(subject_name, sample_name) |> # patient ID, sample ID
  dplyr::mutate(is_duplicate = duplicated(subject_name)) |> # marks 2+ instance of patient ID
  dplyr::mutate(subgroup = dplyr::case_when(stringr::str_detect(sj_associated_diagnoses_disease_code, "G3") ~ "G3",
                                            stringr::str_detect(sj_associated_diagnoses_disease_code, "G4") ~ "G4",
                                            stringr::str_detect(sj_associated_diagnoses_disease_code, "SHH") ~ "SHH",
                                            stringr::str_detect(sj_associated_diagnoses_disease_code, "WNT") ~ "WNT")) |>
  dplyr::rename("sample_accession" = "sample_name") |>
  dplyr::mutate(study = "St. Jude",
                platform = "RNA-seq") |>
  dplyr::select(sample_accession,
                subgroup,
                study,
                is_duplicate,
                platform) |>
  clean_mb_subgroups()

################################################################################
# GSE119926
################################################################################

GSE119926_metadata <- GEOquery::getGEO(filename = GSE119926_metadata_input_filename) |>
  as.data.frame() |>
  dplyr::mutate(is_duplicate = FALSE,
                study = "GSE119926",
                platform =  "Pseudo-bulk",
                is_PDX = dplyr::case_when(stringr::str_detect(source_name_ch1, "patient-derived xenograft") ~ TRUE,
                                          TRUE ~ FALSE)) |>
  dplyr::select(sample_accession = geo_accession,
                title,
                subgroup = methylation.subgroup.ch1,
                study,
                is_duplicate,
                platform,
                is_PDX,
                subtype = methylation.subtype.ch1) |>
  clean_mb_subgroups()

################################################################################
# GSE155446
################################################################################

GSE155446_metadata <- readr::read_csv(GSE155446_metadata_input_filename,
                                      col_types = "c") |>
  dplyr::select(sample_accession = geo_sample_id,
                title = geo_sample_id,
                subgroup = subgroup) |>
  unique() |>
  dplyr::mutate(study = "GSE155446",
                is_duplicate = FALSE,
                platform = "Pseudo-bulk",
                is_PDX = FALSE,
                subtype = NA) |>
  clean_mb_subgroups()

################################################################################
# combine bulk metadata and write to file
################################################################################

# 4 MB subgroups, NA subgroup, Normal subgroup, LGG
dplyr::bind_rows(GSE124814_metadata,
                 GSE164677_metadata,
                 openpbta_mb_metadata,
                 openpbta_lgg_metadata,
                 sj_metadata) |>
  dplyr::filter(!is_duplicate) |>
  readr::write_tsv(file = bulk_metadata_output_filename)

################################################################################
# combine pseudo-bulk metadata and write to file
################################################################################

dplyr::bind_rows(GSE119926_metadata,
                 GSE155446_metadata) |>
  dplyr::filter(!is_duplicate) |>
  readr::write_tsv(file = pseudobulk_metadata_output_filename)
