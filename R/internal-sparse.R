#' Sparse matrix operations
#'
#' @rdname sparse
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param sparse Sparse counts matrix.
#'
#' @return Modified sparse matrix.
#'
#' @description Collapse features with the same feature id by summing them.
#' @param featureids New feature IDs for the collapsed features.
.aggregate_sparse_features <- function(sparse, featureids) {
    rownames(sparse) <- featureids
    sparse <- sparse[!is.na(rownames(sparse)), ]
    aggregate.Matrix(sparse, rownames(sparse), fun = "sum")
}



#' @rdname sparse
#' @description Pool cellular technical replicate counts.
#' @param cellids New cellular ids of the collapsed cells, one for each cell.
.aggregate_sparse_replicates <- function(sparse, cellids) {
    tsparse <- t(sparse)
    rownames(tsparse) <- cellids
    aggregate.Matrix(tsparse, cellids, fun = "sum") %>% t
}



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
