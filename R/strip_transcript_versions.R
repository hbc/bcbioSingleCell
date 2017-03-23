##' Strip transcript versions off of the rownames of a dataframe or matrix
##'
##' @param matrix a matrix or dataframe
##' @return a matrix or dataframe with the transcript versions stripped off
##' @importFrom tools file_path_sans_ext
##' @author Rory Kirchner
##' @export
strip_transcript_versions <- function(matrix) {
    transcripts <- rownames(matrix)
    transcripts <- tools::file_path_sans_ext(transcripts)
    rownames(matrix) <- transcripts
    return(matrix)
}
