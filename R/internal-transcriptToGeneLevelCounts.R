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
    txMap <- tx2gene %>%
        .[.[["enstxp"]] %in% rownames(counts), , drop = FALSE]

    # Detect and handle missing transcript identifiers. These are typically
    # deprecated transcripts in the current Ensembl release, or FASTA
    # spike-in sequences (e.g. EGFP, GAL4). We don't want to simply trash here.
    if (!all(rownames(counts) %in% tx2gene[["enstxp"]])) {
        missing <- rownames(counts) %>%
            .[!. %in% tx2gene[["enstxp"]]] %>%
            sort()
        if (length(missing) <= 200L) {
            # Warn and append unmatched transcripts as genes
            fxn <- warning
            txMap <- data.frame(
                enstxp = missing,
                ensgene = missing,
                row.names = missing,
                stringsAsFactors = FALSE) %>%
                rbind(txMap)
        } else {
            # Stop if there are too many transcript match failures
            fxn <- stop
        }
        fxn(paste(
            length(missing),
            "rows missing from tx2gene:",
            toString(missing)
        ), call. = FALSE)
    }

    message("Converting transcript-level counts to gene-level")
    rownames(counts) <- txMap[["ensgene"]]
    aggregate.Matrix(
        x = counts,
        groupings = rownames(counts),
        fun = "sum")
}
