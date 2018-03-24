#' Transcript to Gene-Level Counts
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @keywords internal
#' @noRd
#'
#' @param object Transcript-level sparse counts matrix (`dgCMatrix`).
#' @param tx2gene Transcript to gene identifier mappings.
#'
#' @return `dgCMatrix`.
.transcriptToGeneLevelCounts <- function(object, tx2gene) {
    assert_is_all_of(object, "dgCMatrix")
    object <- stripTranscriptVersions(object)
    assertIsTx2gene(tx2gene)
    assert_is_subset(rownames(object), rownames(tx2gene))

    # Subset the tx2gene to keep only identifiers present in the matrix
    tx2gene <- tx2gene[rownames(object), , drop = FALSE]
    rownames(tx2gene) <- rownames(object)

    # Detect and handle missing transcript identifiers. These are typically
    # deprecated transcripts in the current Ensembl release, or FASTA
    # spike-in sequences (e.g. EGFP, GAL4). We don't want to simply trash here.
    if (!all(rownames(tx2gene) %in% tx2gene[["txID"]])) {
        match <- tx2gene %>%
            .[!is.na(.[["txID"]]), , drop = FALSE]
        missing <- tx2gene %>%
            .[is.na(.[["txID"]]), , drop = FALSE] %>%
            rownames()
        if (length(missing) > 200L) {
            abort(paste(length(missing), "missing transcripts in tx2gene."))
        }
        warn(paste(
            length(missing), "missing transcripts in tx2gene.",
            "These will be kept but converted to genes:",
            toString(missing)
        ))
        # Warn and append unmatched transcripts as genes
        remap <- data.frame(
            "txID" = missing,
            "geneID" = missing,
            row.names = missing,
            stringsAsFactors = FALSE
        ) %>%
            rbind(match) %>%
            .[rownames(tx2gene), , drop = FALSE]
        tx2gene <- remap
    }

    if (!identical(rownames(object), tx2gene[["txID"]])) {
        abort("Transcript to gene mappings don't match counts matrix rows")
    }

    inform("Converting transcript-level counts to gene-level")
    rownames(object) <- tx2gene[["geneID"]]
    aggregate.Matrix(
        x = object,
        groupings = rownames(object),
        fun = "sum"
    )
}
