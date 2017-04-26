#' Read sample barcodes metadata
#'
#' @author Michael Steinbaugh
#'
#' @param file Sample barcodes metadata file
#'
#' @return Data frame
#' @export
read_sample_barcodes_metadata <- function(file) {
    if (!file.exists(file)) {
        stop("File not found")
    }

    if (grepl("\\.xlsx$", file)) {
        df <- read_excel(file)
    } else {
        df <- read_csv(file)
    }

    # Subset samples
    df <- df %>%
        .[!is.na(.$file_name), ] %>%
        .[!is.na(.$sample_name), ]
    # Reverse complement matches
    df$reverse_complement <- sapply(df$sequence, revcomp)
    # Sample barcode identifier
    df$sample_barcode <- paste(df$file_name,
                               df$reverse_complement,
                               sep = "-")
    return(df)
}
