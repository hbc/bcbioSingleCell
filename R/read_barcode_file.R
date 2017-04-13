#' Load a barcode histogram file generated from bcbio-nextgen
#'
#' @author Rory Kirchner
#'
#' @keywords internal
#'
#' @param filename path to a barcode histogram file
#'
#' @return dataframe of reads per barcode
#' @export
read_barcode_file <- function(filename) {
    read_tsv(filename,
             col_names = c("barcode", "count"),
             progress = FALSE)
}
