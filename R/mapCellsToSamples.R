#' Define Cell to Sample Mappings
#'
#' This function extracts `sampleID` from the `cellID` column using grep
#' matching.
#'
#' @rdname mapCellsToSamples
#' @name mapCellsToSamples
#' @author Michael Steinbaugh
#'
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr bind_rows filter mutate
#' @importFrom magrittr set_colnames
#' @importFrom parallel mclapply
#' @importFrom stringr str_match
#'
#' @param cells Cell identifiers.
#' @param samples Sample identifiers.
#'
#' @return Named factor containing cell IDs as the names, and samples as the
#'   factor levels.
#' @export
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#'
#' cells <- colnames(bcb)
#' samples <- sampleMetadata(bcb)[["sampleID"]]
#'
#' map <- mapCellsToSamples(cells, samples)
#' head(map)
#' levels(map)
mapCellsToSamples <- function(cells, samples) {
    assert_is_character(cells)
    assert_has_no_duplicates(cells)
    assert_is_any_of(samples, c("character", "factor"))
    samples <- unique(as.character(samples))
    list <- bplapply(samples, function(sample) {
        pattern <- paste0("^(", sample, barcodePattern)
        match <- str_match(cells, pattern = pattern) %>%
            as.data.frame() %>%
            set_colnames(c(
                "cellID",
                "sampleID",
                "cellularBarcode",
                # Trailing number is for matching cellranger
                "trailingNumber"
            )) %>%
            select(c("cellID", "sampleID")) %>%
            filter(!is.na(.data[["sampleID"]]))
        vec <- match[["sampleID"]]
        names(vec) <- match[["cellID"]]
        vec
    })
    cell2sample <- unlist(list)
    assert_are_identical(length(cells), length(cell2sample))
    as.factor(cell2sample)
}
