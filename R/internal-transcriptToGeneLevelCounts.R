#' Transcript To Gene-Level Counts
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom dplyr bind_rows
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom magrittr set_rownames
#' @importFrom tibble remove_rownames
#'
#' @param counts Transcript-level sparse counts matrix (`dgCMatrix`).
#' @param tx2gene Transcript to gene identifier mappings.
#'
#' @return `dgCMatrix`.
#' @noRd
.transcriptToGeneLevelCounts <- function(counts, tx2gene) {
    if (!is(counts, "dgCMatrix")) {
        stop("counts matrix must be dgCMatrix class object", call. = FALSE)
    }
    counts <- .stripTranscriptVersions(counts)

    # Subset the tx2gene to keep only identifiers present in the matrix
    map <- tx2gene[rownames(counts), , drop = FALSE]
    rownames(map) <- rownames(counts)

    # Detect and handle missing transcript identifiers. These are typically
    # deprecated transcripts in the current Ensembl release, or FASTA
    # spike-in sequences (e.g. EGFP, GAL4). We don't want to simply trash here.
    if (!all(rownames(map) %in% map[["enstxp"]])) {
        match <- map %>%
            .[!is.na(.[["enstxp"]]), , drop = FALSE]
        missing <- map %>%
            .[is.na(.[["enstxp"]]), , drop = FALSE] %>%
            rownames()
        if (length(missing) > 200L) {
            stop(paste(
                length(missing), "missing transcripts in tx2gene."
            ), call. = FALSE)
        }
        warning(paste(
            length(missing), "missing transcripts in tx2gene.",
            "These will be kept but converted to genes:",
            toString(missing)
        ), call. = FALSE)
        # Warn and append unmatched transcripts as genes
        remap <- data.frame(
            enstxp = missing,
            ensgene = missing,
            row.names = missing,
            stringsAsFactors = FALSE) %>%
            rbind(match) %>%
            .[rownames(map), , drop = FALSE]
        map <- remap
    }

    if (!identical(rownames(counts), map[["enstxp"]])) {
        stop("Transcript to gene mappings don't match counts matrix rows",
             call. = FALSE)
    }

    message("Converting transcript-level counts to gene-level")
    rownames(counts) <- map[["ensgene"]]
    aggregate.Matrix(
        x = counts,
        groupings = rownames(counts),
        fun = "sum")
}
