#' Get transcript-to-gene annotations from
#' \href{http://www.ensembl.org/}{Ensembl}.
#'
#' This function also defines broad classes, which are used in quality control
#' analysis.
#'
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#'
#' @return Data frame.
ensembl_annotations <- function(run) {
    ensembl <- useEnsembl(
        biomart = "ensembl",
        dataset = paste(run$organism, "gene_ensembl", sep = "_"))

    meta <- getBM(
        mart = ensembl,
        attributes = c("ensembl_gene_id",
                       "ensembl_transcript_id",
                       "external_gene_name",
                       "gene_biotype",
                       "chromosome_name")
    ) %>%
        as_tibble %>%
        group_by(!!sym("ensembl_gene_id")) %>%
        arrange(!!sym("ensembl_transcript_id"), .by_group = TRUE)

    # Broad class definitions
    coding <- c("protein_coding")
    decaying <- c("non_stop_decay", "nonsense_mediated_decay")
    noncoding <- c("known_ncrna", "lincRNA", "non_coding")
    srna <- c("miRNA",
              "misc_RNA",
              "ribozyme",
              "rRNA",
              "scaRNA",
              "scRNA",
              "snoRNA",
              "snRNA",
              "sRNA")
    meta$broad_class <-
        case_when(tolower(meta$chromosome_name) == "mt" ~ "mito",
                  # Fix to match Drosophila genome (non-standard)
                  grepl("mito", meta$chromosome_name) ~ "mito",
                  grepl("pseudo", meta$gene_biotype) ~ "pseudo",
                  grepl("TR_", meta$gene_biotype) ~ "TCR",
                  grepl("IG_", meta$gene_biotype) ~ "IG",
                  meta$gene_biotype %in% srna ~ "small",
                  meta$gene_biotype %in% decaying ~ "decaying",
                  meta$gene_biotype %in% noncoding ~ "noncoding",
                  meta$gene_biotype %in% coding ~ "coding",
                  TRUE ~ "other")

    return(meta)
}
