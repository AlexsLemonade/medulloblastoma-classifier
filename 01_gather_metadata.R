# Gather metadata from each data source
#
# Steven Foltz
# November 2022

suppressMessages(library(tidyverse))

data_dir <- here::here("data")
processed_data_dir <- here::here("processed_data")

combined_metadata_output_filename <- file.path(processed_data_dir,
                                               "combined_metadata.tsv")

################################################################################
# functions
################################################################################

clean_mb_subgroups <- function(df){
  
  df <- df %>%
    mutate(subgroup = case_when(subgroup %in% c("E", "Group3", "Group 3", "G3", "Group3_alpha", "Group3_beta", "Group3_gamma", "MB_GRP3") ~ "G3",
                                subgroup %in% c("C", "D", "Group4", "G4", "Group 4", "Group4_alpha", "Group4_beta", "Group4_gamma", "MB_GRP4") ~ "G4",
                                subgroup %in% c("NORM") ~ "Normal",
                                subgroup %in% c("B", "MB_SHH", "SHH", "SHH_alpha", "SHH_beta", "SHH_delta", "SHH_gamma") ~ "SHH",
                                subgroup %in% c("A", "MB_WNT", "WNT", "WNT_alpha", "WNT_beta") ~ "WNT"))

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

# GSE124814 input file names
GSE124814_metadata_input_filename <- file.path(data_dir, "GSE124814",
                                               "GSE124814_sample_descriptions.xlsx")

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
                                        skip = 2) %>%
  separate(title, # separate title into experiment accession and sample id
           into = c("experiment_accession", "id"),
           sep = "-",
           extra = "merge") %>% # split on the first hyphen and keep rest intact
  # filter out E-MTAB-292 per data accessibility
  filter(experiment_accession != "EMTAB292") %>%
  # GSE124814 merged duplicate samples together by averaging expression values
  # and this is indicated in the description column with the word "average"
  mutate(is_duplicate = str_detect(description, "average")) %>%
  # "NA" subgroup corresponds to "Normal"
  mutate(subgroup_supplied_renamed = ifelse(subgroup_supplied_renamed == "NA",
                                            "Normal",
                                            subgroup_supplied_renamed)) %>%
  # treat "Unknown" subgroup as NA -- we need to look further into these
  mutate(subgroup_supplied_renamed = na_if(x = subgroup_supplied_renamed,
                                           y = "Unknown")) %>%
  # isolate the sample accession as its own column ("reanalysis of SAMPLE_ACCESSION (EXPERIMENT_ACCESSION)")
  mutate(sample_accession = word(description, 3)) %>%
  select(sample_accession,
         subgroup_supplied_renamed,
         experiment_accession,
         is_duplicate) %>%
  rename("subgroup" = "subgroup_supplied_renamed",
         "study" = "experiment_accession") %>%
  mutate(platform = "Array") %>%
  clean_mb_subgroups()

################################################################################
# GSE164677
################################################################################

GSE164677_metadata <- read_tsv("data/GSE164677/GSE164677_Asian_MB_RNA-seq.txt.gz",
                               col_names = FALSE,
                               col_types = "c",
                               n_max = 2) %>%
  column_to_rownames(var = "X1") %>%
  t() %>%
  as_tibble() %>%
  rename("sample_accession" = "gene",
         "subgroup" = "group") %>%
  select(sample_accession,
         subgroup) %>%
  mutate(study = "GSE164677",
         is_duplicate = FALSE,
         platform = "RNA-seq") %>%
  clean_mb_subgroups()

################################################################################
# OpenPBTA (MB)
################################################################################

openpbta_mb_metadata <- read_tsv(file = "data/OpenPBTA/pbta-histologies.tsv",
                                 col_types = "c") %>%
  filter(experimental_strategy == "RNA-Seq",
         short_histology == "Medulloblastoma") %>%
  separate(molecular_subtype, # separate molecular_subtype into molecular and subgroup
           into = c("molecular", "subgroup"),
           sep = ", ",
           extra = "merge") %>%
  mutate(subgroup = na_if(x = subgroup,
                          y = "To be classified")) %>% 
  arrange(Kids_First_Participant_ID) %>% # patient ID
  mutate(is_duplicate = duplicated(Kids_First_Participant_ID)) %>% # marks 2+ instance of patient ID
  rename("sample_accession" = "Kids_First_Biospecimen_ID") %>%
  mutate(study = "OpenPBTA",
         platform = "RNA-seq") %>%
  select(sample_accession,
         subgroup,
         study,
         is_duplicate,
         platform) %>%
  clean_mb_subgroups()

################################################################################
# OpenPBTA (LGG)
################################################################################

openpbta_lgg_metadata <- read_tsv(file = "data/OpenPBTA/pbta-histologies.tsv",
                                  col_types = "c") %>%
  filter(experimental_strategy == "RNA-Seq",
         pathology_diagnosis == "Low-grade glioma/astrocytoma (WHO grade I/II)",
         short_histology == "LGAT") %>%
  mutate(subgroup = "LGG") %>%
  arrange(Kids_First_Participant_ID) %>% # patient ID
  mutate(is_duplicate = duplicated(Kids_First_Participant_ID)) %>% # marks 2+ instance of patient ID
  rename("sample_accession" = "Kids_First_Biospecimen_ID") %>%
  mutate(study = "OpenPBTA",
         platform = "RNA-seq") %>%
  select(sample_accession,
         subgroup,
         study,
         is_duplicate,
         platform) # do not clean_mb_subgroups()

################################################################################
# St. Jude
################################################################################

sj_metadata <- read_tsv("data/stjudecloud/SAMPLE_INFO.txt",
                        col_types = "c") %>%
  separate(sj_embargo_date,
           into = c("embargo_mon",
                    "embargo_day",
                    "embargo_year"),
           sep = "-") %>%
  mutate(embargo_year = as.numeric(embargo_year)) %>%
  filter(embargo_year < 2023) %>% # keep samples with embargo ending before 2023
  arrange(subject_name) %>% # patient ID
  mutate(is_duplicate = duplicated(subject_name)) %>% # marks 2+ instance of patient ID
  mutate(subgroup = case_when(str_detect(sj_associated_diagnoses_disease_code, "G3") ~ "G3",
                              str_detect(sj_associated_diagnoses_disease_code, "G4") ~ "G4",
                              str_detect(sj_associated_diagnoses_disease_code, "SHH") ~ "SHH",
                              str_detect(sj_associated_diagnoses_disease_code, "WNT") ~ "WNT")) %>%
  rename("sample_accession" = "sample_name") %>%
  mutate(study = "St. Jude",
         platform = "RNA-seq") %>%
  select(sample_accession,
         subgroup,
         study,
         is_duplicate,
         platform) %>%
  clean_mb_subgroups()

################################################################################
# combine metadata and write to file
################################################################################

# 4 MB subgroups, NA subgroup, Normal subgroup, LGG
bind_rows(GSE124814_metadata,
          GSE164677_metadata,
          openpbta_mb_metadata,
          openpbta_lgg_metadata,
          sj_metadata) %>%
  filter(!is_duplicate) %>% 
  write_tsv(file = combined_metadata_output_filename)
