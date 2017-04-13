#' Load a CSV file with readr, setting a given column as the rownames
#'
#' @author Rory Kirchner
#'
#' @keywords internal
#'
#' @param filename CSV to read
#' @param column column to make into the rownames
#'
#' @return Data frame
#' @export
load_csv_with_rownames <- function(filename, column) {
    dat <- read_csv(filename, progress = FALSE) %>%
        as.data.frame
    rownames(dat) <- dat[, column]
    dat[, column] <- NULL
    return(dat)
}
