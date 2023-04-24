# Downloading gene expression feature counts from St. Jude Cloud

[St. Jude Cloud](https://www.stjude.cloud/) requires users to create an account before downloading data from the [Data Browser](https://platform.stjude.cloud/data/diseases/paired-tumor-normal).
Those who are not St. Jude employees need to first create a DNAnexus account since data are hosted by DNAnexus.

Once registered, login and use the following "Share Selection" URL to locate RNA-seq feature counts data from medulloblastoma (MBL) diagnosis samples:

`https://platform.stjude.cloud/data/diseases/tumor?file_type=FEATURE_COUNTS&seq_type=RNA-SEQ&sample_type=diagnosis&search=MBL&selected_tags=AMBLNWSG3,AMBLNWSG4,DMBLSHH,AMBL,AMBLSHH,MBL,MBLG3,MBLG4,MBLSHH,MBLWNT`

A set of 123 samples downloaded for this study (accessed 2022-06-16) is available at [data/stjudecloud_sample_names.tsv](data/stjudecloud_sample_names.tsv), with 57 samples matching the following criteria:

- Embargo date before 2023
- Not a duplicate sample (drop subsequent samples from the same patient)
- Has a well-defined medulloblastoma subgroup
