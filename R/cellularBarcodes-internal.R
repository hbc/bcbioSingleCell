#' Proportional Cellular Barcodes
#'
#' @rdname proportionalCB-internal
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @details Modified version of Klein Lab MATLAB code.
#'
#' @return [tibble].
.proportionalCB <- function(object) {
    metadata <- metadata(object)[["sampleMetadata"]] %>%
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



#' Read Cellular Barcode File
#'
#' @rdname readCBFile-internal
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param file Cellular barcode TSV file.
#'
#' @return [tibble].
.readCBFile <- function(file) {
    readFileByExtension(
        file,
        col_names = c("cellularBarcode", "reads"),
        col_types = "ci") %>%
        mutate(cellularBarcode = str_replace_all(
            .data[["cellularBarcode"]], "-", "_"))
}



#' Cellular Barcodes List
#'
#' @rdname cbList-internal
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param sampleDirs Sample directories.
#'
#' @return [list].
.cbList <- function(sampleDirs) {
    files <- sampleDirs %>%
        file.path(paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        set_names(names(sampleDirs))
    if (!all(file.exists(files))) {
        stop("Cellular barcode file missing")
    }
    message("Reading cellular barcode distributions")
    pblapply(seq_along(files), function(a) {
        .readCBFile(files[a])
    }) %>% set_names(names(sampleDirs))
}



#' Bind Cellular Barcodes
#'
#' @rdname bindCB-internal
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param list List of cellular barcodes.
#'
#' @return [tibble].
.bindCB <- function(list) {
    lapply(seq_along(list), function(a) {
        sampleID <- names(list)[[a]] %>% camel
        list[[a]] %>%
            mutate(sampleID = !!sampleID)
    }) %>%
        as("tibble") %>%
        bind_rows %>%
        mutate(rowname = paste(.data[["sampleID"]],
                               .data[["cellularBarcode"]],
                               sep = "_"))
}
