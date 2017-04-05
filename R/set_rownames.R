#' Set row names
#'
#' Wrapper function for \code{rownames} that works in a chain operation
#'
#' @author Michael Steinbaugh
#'
#' @import tibble
#'
#' @param df Data frame
#' @param column Column to use for \code{rownames}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set_rownames(df, "description")
#' }
set_rownames <- function(df, column) {
    # Setting rownames on a tibble is deprecated. Therefore, we must coerce to
    # data frame first, if necessary.
    if (tibble::is_tibble(df)) {
        df <- as.data.frame(df)
    }
    rownames(df) <- df[[column]]
    return(df)
}
