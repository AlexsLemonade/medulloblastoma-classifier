# Functions to convert gene expression counts to TPM (transcripts per million) values
#
# Steven Foltz
# February 2023
#
# Sources:
# 
# St. Jude RNA-seq workflow v2.0.0 https://stjudecloud.github.io/rfcs/0001-rnaseq-workflow-v2.0.0.html
# which uses htseq-count with https://www.gencodegenes.org/human/release_31.html to create gene level counts.
# One approach to deriving gene lengths: https://www.biostars.org/p/83901/
# Description of TPM values: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

get_GENCODE_gene_lengths <- function(gtf_filepath, GENCODE_gene_lengths_filepath) {
  # If a GENCODE gene length file does not already exist, create one using GTF
  # 
  # Inputs
  #   gtf_filepath: file path to GTF file (can be .gz)
  #   GENCODE_gene_lengths_filepath: GENCODE gene length file path to read from or write to
  #
  # Outputs
  #   Returns a GENCODE gene length data frame
  #   May write GENCODE gene length data frame to GENCODE_gene_lengths_filepath
  
  if (file.exists(GENCODE_gene_lengths_filepath) ) {
    
    message(glue::glue(c(GENCODE_gene_lengths_filepath),
                       " already exists and will not be re-created."))
    
    GENCODE_gene_lengths_df <- readr::read_tsv(GENCODE_gene_lengths_filepath,
                                              show_col_types = FALSE)
    
  } else {
    
    # follows example from https://www.biostars.org/p/83901/
    
    # create transcriptome database from GTF file
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_filepath, format = "gtf")
    
    # gather exons from each transcript at gene level and
    # collapse each gene's exons to form union of bases covered
    exonic.gene.sizes <- GenomicFeatures::exonsBy(txdb, by = "gene") %>%
      GenomicRanges::reduce() %>%
      IRanges::width() %>%
      sum()
    
    # create data frame of gene lengths
    # ENSEMBL IDs are given with a version, like .1 or .2,
    #   but we need to remove the version for mapping to the gene expression data frame
    # Some ENSEMBLE IDs have a version like .20_PAR_Y,
    #   which all have duplicate ENSGs that have regular versions like .20,
    #   so we can drop PAR_Y versions and keep regular versions to avoid multi-mapping ENSGs
    GENCODE_gene_lengths_df <- dplyr::tibble(ENSEMBL_with_version = names(exonic.gene.sizes),
                                             gene_length = as.vector(exonic.gene.sizes)) %>%
      tidyr::separate(ENSEMBL_with_version,
                      into = c("ENSEMBL", "version"),
                      sep = "\\.") %>%
      dplyr::filter(!stringr::str_detect(version, "_PAR_Y")) %>% # remove Y paralogs
      select(ENSEMBL, gene_length) %>%
      readr::write_tsv(file = GENCODE_gene_lengths_filepath)    
    
  }
  
  return(GENCODE_gene_lengths_df)

}

convert_gene_counts_to_TPM <- function(genex_df, gene_lengths_df){
  
  # Given a data frame of gene counts and the length of each gene, convert to TPM
  #
  # Inputs
  #   genex_df: gene expression data frame with gene names as row names
  #   gene_lengths_df: gene length data frame with columns ENSEMBL and gene_length
  #
  # Outputs
  #   Returns a gene expression matrix with TPM values
  
  # reduce gene expression df to those genes with a known length
  genex_df <- genex_df[rownames(genex_df) %in% gene_lengths_df$ENSEMBL,]
  
  # harmonize the order of genes and pull out a vector gene lengths
  gene_length_vector <- tibble::tibble(genex_genes = rownames(genex_df)) %>%
    dplyr::left_join(gene_lengths_df,
                     by = c("genex_genes" = "ENSEMBL")) %>%
    pull(gene_length)
  
  # TPM = reads_per_kilobase/scaling_factor
  # where reads_per_kilobase is the read count divided by length of gene in kilobases
  # and scaling factor makes each column sum to 1e6 (sum the reads per KB and divide by 1e6)
  # (there is one scaling factor per column)
  
  # calculate reads per kilobase
  reads_per_kilobase_df <- apply(genex_df, 2, function(x) 1000*x/gene_length_vector)
  
  # calculate scaling factor
  per_million_scaling_factor <- apply(reads_per_kilobase_df, 2, function(x) sum(x)/1e6)
  
  # divide the values in each row by the corresponding scaling factor
  # transpose back to get genes x samples shape
  tpm_df <- as.data.frame(t(apply(reads_per_kilobase_df, 1, function(x) x/per_million_scaling_factor)))
  
  return(tpm_df)
  
}
