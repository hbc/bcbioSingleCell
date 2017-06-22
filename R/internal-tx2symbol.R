# TODO Switch to annotable method instead of using fromto

#' Convert a transcript matrix to a gene symbol matrix
#'
#' @rdname tx2symbol
#' @keywords internal
#'
#' @author Rory Kirchner
#'
#' @param matrix Counts matrix.
#' @param strip Strip transcript versions from the matrix before aggregating.
#'
#' @return Matrix of counts aggregated by summing the "to" ids.
.tx2symbol <- function(matrix, strip = TRUE) {
    # TODO Add grep match on identifiers, to see if stripping is necessary.
    # Then we can remove this user-defined parameter.
    if (isTRUE(strip)) {
        matrix <- strip_transcript_versions(matrix)
    }
    lookup <- data.frame(from = rownames(matrix)) %>%
        left_join(fromto, by = "from")
    rownames(matrix) <- lookup$to
    matrix <- matrix[!is.na(rownames(matrix)), ]
    aggregate.Matrix(matrix, row.names(matrix), fun = "sum")
}
