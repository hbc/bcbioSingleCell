#' Read a CSV file with readr, setting a given column as the rownames.
#'
#' @keywords internal
#' @author Rory Kirchner
#'
#' @param filename CSV to read.
#' @param column Column to make into the rownames.
#'
#' @return Data frame.
#' @export
read_csv_with_rownames <- function(filename, column) {
    dat <- read_csv(filename, progress = FALSE) %>%
        as.data.frame
    rownames(dat) <- dat[, column]
    dat[, column] <- NULL
    return(dat)
}
