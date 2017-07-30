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
    cb <- bcbio(object, "cellularBarcodes")
    lapply(seq_along(cb), function(a) {
        cb <- cb[[a]] %>%
            mutate(log10Reads = log10(.data[["reads"]]))
        cbHist <- hist(cb[["log10Reads"]], n = 100L, plot = FALSE)
        # `fLog` in Klein Lab code
        counts <- cbHist[["counts"]]
        # `xLog` in Klein Lab code
        mids <-  cbHist[["mids"]]
        tibble(
            sampleID = names(cb)[[a]],
            log10ReadsPerCell = mids,
            proportionOfCells = counts * (10L ^ mids) /
                sum(counts * (10L ^ mids)))
    }) %>%
        set_names(names(cb)) %>%
        bind_rows %>%
        left_join(metadata, by = "sampleID")
}
