#' @rdname sparse
#' @description Strip transcript versions
#' @return Sparse counts matrix with the transcript versions stripped.
.strip_transcript_versions <- function(sparse) {
    check_sparse(sparse)
    transcripts <- rownames(sparse)
    # Tight pattern matching against Ensembl stable transcript IDs
    # http://www.ensembl.org/info/genome/stable_ids/index.html
    pattern <- "^(ENS[A-Z]{3}T\\d{11})\\.\\d+$"
    if (any(str_detect(transcripts, pattern))) {
        transcripts <- str_replace(transcripts, pattern, "\\1")
    }
    rownames(sparse) <- transcripts
    sparse
}
