#' Read sample barcodes metadata
#'
#' @author Michael Steinbaugh
#'
#' @param file Sample barcodes metadata file
#'
#' @return Data frame
#' @export
read_metadata <- function(file) {
    if (!file.exists(file)) {
        stop("File not found")
    }

    if (grepl("\\.xlsx$", file)) {
        metadata <- read_excel(file)
    } else {
        metadata <- read_csv(file)
    }

    # Subset samples
    metadata <- metadata %>%
        .[!is.na(.$file_name), ] %>%
        .[!is.na(.$sample_name), ]
    # Reverse complement matches
    metadata$reverse_complement <- sapply(metadata$sequence, revcomp)
    # Sample barcode identifier
    metadata$sample_barcode <- paste(metadata$file_name,
                                     metadata$reverse_complement,
                                     sep = "-")
    return(metadata)
}
