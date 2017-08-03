#' Proportional Cellular Barcodes
#'
#' @rdname proportionalCB-internal
#' @family Cellular Barcode Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @details Modified version of Klein Lab MATLAB code.
#'
#' @return [tibble].
.proportionalCB <- function(object) {
    metadata <- sampleMetadata(object) %>%
        .[, metaPriorityCols]
    lst <- bcbio(object, "cellularBarcodes")
    lapply(seq_along(lst), function(a) {
        cb <- lst[[a]] %>%
            mutate(log10Reads = log10(.data[["reads"]]))
        cbHist <- hist(cb[["log10Reads"]], n = 100L, plot = FALSE)
        # `fLog` in Klein Lab code
        counts <- cbHist[["counts"]]
        # `xLog` in Klein Lab code
        mids <-  cbHist[["mids"]]
        tibble(
            sampleID = names(lst)[[a]],
            # log10 read counts per cell
            log10Count = mids,
            # Proportion of cells
            proportion = counts * (10L ^ mids) /
                sum(counts * (10L ^ mids)))
    }) %>%
        set_names(names(lst)) %>%
        bind_rows %>%
        left_join(metadata, by = "sampleID")
}
