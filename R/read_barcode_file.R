##' Load a barcode histogram file generated from bcbio-nextgen
##'
##' @param filename path to a barcode histogram file
##' @return dataframe of reads per barcode
##' @importFrom readr read_tsv
##' @keywords internal
##' @author Rory Kirchner
##' @export
read_barcode_file <- function(filename) {
    readr::read_tsv(filename,
                    col_names = c("barcode", "count"),
                    progress = FALSE) %>%
        return
}
