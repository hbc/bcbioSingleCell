#' Strip Transcript Versions
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param sparseCounts Sparse counts matrix (`dgCMatrix`).
#'
#' @return `dgCMatrix`.
.stripTranscriptVersions <- function(sparseCounts) {
    transcripts <- rownames(sparseCounts)
    # Pattern matching against Ensembl transcript IDs
    # http://www.ensembl.org/info/genome/stable_ids/index.html
    # Examples: ENST (human); ENSMUST (mouse)
    enstxpPattern <- "^(ENS.*T\\d{11})\\.\\d+$"
    if (any(grepl(x = transcripts, pattern = enstxpPattern))) {
        transcripts <- gsub(
            x = transcripts,
            pattern = enstxpPattern,
            replacement = "\\1")
        rownames(sparseCounts) <- transcripts
    }
    if (any(grepl(x = transcripts, pattern = "\\.\\d+$"))) {
        abort("Incomplete transcript version removal")
    }
    sparseCounts
}
