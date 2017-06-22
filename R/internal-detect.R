#' Detect scRNA-seq analysis pipeline
#'
#' Currently supports bcbio-nextgen and 10X Chromium.
#'
#' @rdname detect_pipeline
#'
#' @param upload_dir Final upload directory.
#'
#' @author Michael Steinbaugh



.detect_tenx <- function(upload_dir) {
    # Recursively find the counts matrix file
    sample_dirs <- list.files(
        upload_dir,
        full.names = TRUE,
        pattern = "matrix.mtx",
        recursive = TRUE) %>%
        normalizePath %>%
        dirname %>%
        set_names(basename(.))
    if (!length(sample_dirs)) {
        return(FALSE)
    }

    # Now check for the presence of required TSV files
    barcodes_files <- file.path(sample_dirs, "barcodes.tsv")
    if (!all(file.exists(barcodes_files))) {
        return(FALSE)
    }

    genes_files <- file.path(sample_dirs, "genes.tsv")
    if (!all(file.exists(genes_files))) {
        return(FALSE)
    }

    TRUE
}
