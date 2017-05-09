#' Load a barcode histogram file generated from bcbio-nextgen
#'
#' @author Rory Kirchner
#' @keywords internal
#'
#' @param file_name path to a barcode histogram file
#'
#' @return dataframe of reads per barcode
#' @export
read_barcode_file <- function(file_name) {
    read_tsv(file_name,
             col_names = c("barcode", "count"),
             progress = FALSE)
}
