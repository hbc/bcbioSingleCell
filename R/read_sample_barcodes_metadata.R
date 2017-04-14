#' Read sample barcodes metadata
#'
#' @author Michael Steinbaugh
#'
#' @param file Sample barcodes metadata file
#'
#' @return Data frame
#' @export
read_sample_barcodes_metadata <- function(
    file = file.path("meta", "sample_barcodes.xlsx")) {
    if (file.exists(file)) {
        if (grepl("\\.xlsx$", file)) {
            df <- read_excel(file)
        } else {
            df <- read_csv(file)
        }
        # Reverse complement matches
        df$reverse_complement <- sapply(df$sequence, revcomp)
        # Unique identifier for each barcode
        df$sample_barcode <- paste(df$samplename,
                                   df$reverse_complement,
                                   sep = "-")
        return(df)
    } else {
        return(NULL)
    }
}
