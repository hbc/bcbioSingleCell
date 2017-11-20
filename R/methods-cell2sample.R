#' Cell to Sample Mappings
#'
#' @rdname cell2sample
#' @name cell2sample
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
NULL



# Constructors ====
#' Define Cell to Sample Mappings
#'
#' This function extracts `sampleID` from the `cellID` column using grep
#' matching.
#'
#' @importFrom dplyr bind_rows filter mutate
#' @importFrom magrittr set_colnames
#' @importFrom parallel mclapply
#' @importFrom stringr str_match
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
.cell2sample <- function(cells, samples) {
    cells <- as.character(cells)
    samples <- as.character(samples)
    list <- mclapply(seq_along(samples), function(a) {
        pattern <- paste0("^(", samples[[a]], barcodePattern)
        match <- str_match(cells, pattern = pattern) %>%
            as.data.frame() %>%
            set_colnames(
                c("cellID",
                  "sampleID",
                  "cellularBarcode",
                  # Trailing number is for matching cellranger
                  "trailingNumber")
            ) %>%
            select(c("cellID", "sampleID")) %>%
            filter(!is.na(.data[["sampleID"]]))
        match
    })
    cell2sample <- bind_rows(list) %>%
        mutate(sampleID = as.factor(.data[["sampleID"]]))
    if (!identical(length(cells), nrow(cell2sample))) {
        stop("Failed to correctly match sample IDs to cellular barcodes",
             call. = FALSE)
    }
    cell2sample
}



# Methods ====
#' @rdname cell2sample
#' @param samples Sample identifiers.
#' @export
setMethod(
    "cell2sample",
    signature("character"),
    function(object, samples) {
        .cell2sample(
            cells = object,
            samples = samples
        )
    })



#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("bcbioSingleCell"),
    function(object) {
        # Define the cell2sample mappings
        # This uses a stashed `data.frame` as of v0.0.22, for better speed
        cell2sample <- metadata(object)[["cell2sample"]]
        if (is.null(cell2sample)) {
            cell2sample <- .cell2sample(
                cells = colData[["cellID"]],
                samples = metadata[["sampleID"]]
            )
        }
        cell2sample[["sampleID"]] <- as.factor(cell2sample[["sampleID"]])
        cell2sample
    })
