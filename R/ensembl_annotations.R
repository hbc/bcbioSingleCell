#' Get transcript annotations from Ensembl
#'
#' This function also defines broad classes, which are used in quality control
#' analysis
#'
#' @author Michael Steinbaugh
#' @keywords metadata
#'
#' @import dplyr
#' @importFrom biomaRt getBM useMart
#'
#' @param organism Organism identifier
#'
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' ensembl_annotations("mmusculus")
#' }
ensembl_annotations <- function(organism) {
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                             paste(organism, "gene_ensembl", sep = "_"))
    # mart_attributes <- biomaRt::listAttributes(mart)

    df <- biomaRt::getBM(mart = mart,
                         attributes = c("ensembl_transcript_id",
                                        "ensembl_gene_id",
                                        "external_transcript_name",
                                        "external_gene_name",
                                        "transcript_biotype",
                                        "gene_biotype",
                                        "chromosome_name")) %>%
        dplyr::arrange_(.dots = "ensembl_transcript_id")
    rownames(df) <- df[[1]]

    # Broad class definitions
    coding <- c("protein_coding")
    decaying <- c("non_stop_decay",
                  "nonsense_mediated_decay")
    noncoding <- c("known_ncrna",
                   "lincRNA",
                   "non_coding")
    srna <- c("miRNA",
              "misc_RNA",
              "ribozyme",
              "rRNA",
              "scaRNA",
              "scRNA",
              "snoRNA",
              "snRNA",
              "sRNA")
    df$broad_class <-
        dplyr::case_when(grepl("mito", df$chromosome_name) ~ "mito",
                         grepl("pseudo", df$gene_biotype) ~ "pseudo",
                         grepl("TR_", df$gene_biotype) ~ "TCR",
                         grepl("IG_", df$gene_biotype) ~ "IG",
                         df$gene_biotype %in% srna ~ "small",
                         df$gene_biotype %in% decaying ~ "decaying",
                         df$gene_biotype %in% noncoding ~ "noncoding",
                         df$gene_biotype %in% coding ~ "coding",
                         TRUE ~ "other")

    return(df)
}
