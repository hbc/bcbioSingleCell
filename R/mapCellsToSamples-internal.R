#' Define Cell to Sample Mappings
#'
#' This function extracts `sampleID` from the `cellID` column using grep
#' matching.
#'
#' @name mapCellsToSamples
#' @family Data Functions
#' @author Michael Steinbaugh
#' @noRd
#'
#' @param cells `character`. Cell identifiers.
#' @param samples `character`. Sample identifiers.
#'
#' @return Named `factor` containing cells as the names, and samples as the
#'   factor levels.
#'
#' @examples
#' # bcbioSingleCell ====
#' object <- indrops_small
#'
#' # Cell IDs are the column names of object
#' cells <- colnames(object)
#' print(cells)
#'
#' # Use the sample ID for the mappings, which must be present in the cell ID
#' samples <- names(sampleNames(object))
#' # Alternate approach:
#' # samples <- rownames(sampleData(object))
#' # Sample ID must be the prefix of the cell IDs
#' stopifnot(all(grepl(paste0("^", samples[[1L]]), cells)))
#' print(samples)
#'
#' x <- .mapCellsToSamples(cells = cells, samples = samples)
#' glimpse(x)
.mapCellsToSamples <- function(cells, samples) {
    assert_is_character(cells)
    assert_has_no_duplicates(cells)
    assert_is_any_of(samples, c("character", "factor"))
    samples <- unique(as.character(samples))

    # Early return if `cells` don't have a separator and `samples` is a string.
    # Shouldn't happen normally but is necessary for some Seurat datasets.
    if (!any(grepl("[_-]", cells)) && is_a_string(samples)) {
        cell2sample <- factor(replicate(n = length(cells), expr = samples))
        names(cell2sample) <- cells
        return(cell2sample)
    }

    list <- lapply(samples, function(sample) {
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
            .[, c("cellID", "sampleID")] %>%
            .[!is.na(.[["sampleID"]]), , drop = FALSE]
        vec <- match[["sampleID"]]
        names(vec) <- match[["cellID"]]
        vec
    })

    cell2sample <- unlist(list)
    assert_are_identical(length(cells), length(cell2sample))
    as.factor(cell2sample)
}
