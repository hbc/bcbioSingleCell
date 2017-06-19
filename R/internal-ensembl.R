#' Ensembl annotations
#'
#' Download transcript-to-gene annotations from
#' [Ensembl](http://www.ensembl.org/). This function also defines broad classes,
#' which are used in quality control analysis.
#'
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#'
#' @return Data frame.
ensembl_annotations <- function(run) {
    # Broad class definitions
    coding <- c("protein_coding")
    decaying <- c("non_stop_decay", "nonsense_mediated_decay")
    noncoding <- c("known_ncrna", "lincRNA", "non_coding")
    srna <- c("miRNA", "misc_RNA", "ribozyme", "rRNA", "scaRNA", "scRNA",
              "snoRNA", "snRNA", "sRNA")
    ensembl <- useEnsembl(
        biomart = "ensembl",
        dataset = paste(run$organism, "gene_ensembl", sep = "_"))
    getBM(mart = ensembl,
          attributes = c("ensembl_gene_id",
                         "ensembl_transcript_id",
                         "external_gene_name",
                         "gene_biotype",
                         "chromosome_name")) %>%
        as_tibble %>%
        group_by(!!sym("ensembl_gene_id")) %>%
        arrange(!!sym("ensembl_transcript_id"), .by_group = TRUE) %>%
        mutate(broad_class = case_when(
            tolower(.data$chromosome_name) == "mt" ~ "mito",
            # Fix to match Drosophila genome (non-standard)
            grepl("mito", .data$chromosome_name) ~ "mito",
            grepl("pseudo", .data$gene_biotype) ~ "pseudo",
            grepl("TR_", .data$gene_biotype) ~ "TCR",
            grepl("IG_", .data$gene_biotype) ~ "IG",
            .data$gene_biotype %in% srna ~ "small",
            .data$gene_biotype %in% decaying ~ "decaying",
            .data$gene_biotype %in% noncoding ~ "noncoding",
            .data$gene_biotype %in% coding ~ "coding",
            TRUE ~ "other"))
}



#' Metadata table for knit report.
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#'
#' @return \code{\link[knitr]{kable}}.
#' @export
metadata_table <- function(run) {
    run$metadata %>%
        as_tibble %>%
        remove_rownames %>%
        mutate(file_name = NULL) %>%
        set_names(str_replace_all(colnames(.), "_", " ")) %>%
        kable(caption = "Sample metadata")
}



#' Convert a transcript matrix to a gene symbol matrix
#'
#' @keywords internal
#' @author Rory Kirchner
#'
#' @param matrix Counts matrix.
#' @param fromto Data frame where "from" is the original id and "to" is the id
#'   to convert to.
#' @param strip Strip transcript versions from the matrix before aggregating.
#'
#' @return Matrix of counts aggregated by summing the "to" ids.
#' @export
tx2symbol <- function(matrix, fromto, strip = FALSE) {
    if (strip) {
        matrix <- strip_transcript_versions(matrix)
    }
    lookup <- data.frame(from = rownames(matrix)) %>%
        left_join(fromto, by = "from")
    rownames(matrix) <- lookup$to
    matrix <- matrix[!is.na(rownames(matrix)), ]
    aggregate.Matrix(matrix, row.names(matrix), fun = "sum")
}
