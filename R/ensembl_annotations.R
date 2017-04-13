#' Get transcript annotations from Ensembl
#'
#' This function also defines broad classes, which are used in quality control
#' analysis
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run
#'
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' ensembl_annotations(run)
#' }
ensembl_annotations <- function(run) {
    check_run(run)

    # version = "87"
    ensembl <- useEnsembl(
        biomart = "ensembl",
        dataset = paste(run$organism, "gene_ensembl", sep = "_"),
        host = "dec2016.archive.ensembl.org")

    df <- getBM(
        mart = ensembl,
        attributes = c("ensembl_transcript_id",
                       "external_gene_name",
                       "gene_biotype",
                       "chromosome_name")
    ) %>%
        arrange_(.dots = "ensembl_transcript_id") %>%
        set_rownames(.$ensembl_transcript_id)

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
        case_when(tolower(df$chromosome_name) == "mt" ~ "mito",
                  # Fix to match Drosophila genome (non-standard)
                  grepl("mito", df$chromosome_name) ~ "mito",
                  grepl("pseudo", df$gene_biotype) ~ "pseudo",
                  grepl("TR_", df$gene_biotype) ~ "TCR",
                  grepl("IG_", df$gene_biotype) ~ "IG",
                  df$gene_biotype %in% srna ~ "small",
                  df$gene_biotype %in% decaying ~ "decaying",
                  df$gene_biotype %in% noncoding ~ "noncoding",
                  df$gene_biotype %in% coding ~ "coding",
                  TRUE ~ "other")
    df$chromosome_name <- NULL

    return(df)
}
