#' Load a barcode histogram file
#'
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param file_name Barcode histogram TSV file.
#'
#' @return Named numeric vector.
#' @export
read_barcode_file <- function(file_name) {
    df <- read_tsv(file_name,
                   col_names = c("cellular_barcode", "reads"),
                   progress = FALSE)
    df$reads %>%
        as.numeric %>%
        set_names(df$cellular_barcode) %>%
        sort(decreasing = TRUE)
}
