#' Strip transcript versions off of the rownames of a dataframe or matrix
#'
#' @author Rory Kirchner
#'
#' @keywords internal
#'
#' @param matrix a matrix or dataframe
#'
#' @return a matrix or dataframe with the transcript versions stripped off
#' @export
strip_transcript_versions <- function(matrix) {
    transcripts <- rownames(matrix)
    transcripts <- file_path_sans_ext(transcripts)
    rownames(matrix) <- transcripts
    return(matrix)
}
