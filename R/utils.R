##' Load a CSV file with readr, setting a given column as the rownames
##'
##' @param filename CSV to read
##' @param column column to make into the rownames
##' @return dataframe
##' @importFrom readr read_csv
##' @author Rory Kirchner
load_csv_with_rownames <- function(filename, column) {
    dat <- read_csv(filename, progress=FALSE) %>%
        as.data.frame()
    rownames(dat) <- dat[, column]
    dat[, column] <- NULL
    return(dat)
}
