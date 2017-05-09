#' Strip transcript versions
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param counts Sparse counts matrix
#'
#' @return Sparse counts matrix with the transcript versions stripped
#' @export
strip_transcript_versions <- function(counts) {
    check_sparse(counts)
    transcripts <- rownames(counts)
    # Tight pattern matching against Ensembl stable transcript IDs
    # http://www.ensembl.org/info/genome/stable_ids/index.html
    pattern <- "^(ENS[A-Z]{3}T\\d{11})\\.\\d+$"
    if (any(str_detect(transcripts, pattern))) {
        transcripts <- str_replace(transcripts, pattern, "\\1")
    }
    rownames(counts) <- transcripts
    return(counts)
}
