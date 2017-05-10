#' Convert a transcript matrix to a gene symbol matrix
#'
#' @keywords internal
#' @author Rory Kirchner
#'
#' @param matrix a matrix of counts
#' @param fromto a dataframe where "from" is the original id and "to" is the id
#'   to convert to
#' @param strip strip transcript versions from the matrix before aggregating
#'
#' @return a matrix of counts aggregated by summing the "to" ids
#' @export
tx2symbol <- function(matrix, fromto, strip = FALSE) {
    if (strip) {
        matrix <- strip_transcript_versions(matrix)
    }
    lookup <- data.frame(from = rownames(matrix)) %>%
        left_join(fromto, by = "from")
    rownames(matrix) <- lookup$to
    matrix <- matrix[!is.na(rownames(matrix)), ]
    matrix <- aggregate.Matrix(matrix,
                               row.names(matrix),
                               fun = "sum")
    return(matrix)
}
