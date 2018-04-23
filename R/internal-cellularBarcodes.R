#' Bind Cellular Barcodes
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param list List of cellular barcodes.
#'
#' @return `data.frame`.
.bindCellularBarcodes <- function(list) {
    assert_is_list(list)
    mcmapply(
        sampleID = names(list),
        x = list,
        FUN = function(x, sampleID) {
            assert_is_data.frame(x)
            # Add the sampleID as a column
            x[["sampleID"]] <- sampleID
            x
        },
        USE.NAMES = TRUE,
        SIMPLIFY = FALSE
    ) %>%
        bind_rows() %>%
        mutate(
            rowname = paste(
                !!sym("sampleID"),
                !!sym("cellularBarcode"),
                sep = "_"
            ),
            sampleID = as.factor(!!sym("sampleID"))
        ) %>%
        as.data.frame() %>%
        column_to_rownames() %>%
        # Reorder the columns before return
        select(!!!syms(c("sampleID", "cellularBarcode", "nCount"))
}



#' Cellular Barcodes List
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param sampleDirs Sample directories.
#'
#' @return `list`.
.cellularBarcodesList <- function(sampleDirs) {
    files <- file.path(
        normalizePath(sampleDirs, winslash = "/", mustWork = TRUE),
        paste(basename(sampleDirs), "barcodes.tsv", sep = "-")
    )
    names(files) <- names(sampleDirs)
    assert_all_are_existing_files(files)
    message("Reading cellular barcode distributions")
    list <- mclapply(files, function(file) {
        read_tsv(
            file = file,
            col_names = c("cellularBarcode", "nCount"),
            col_types = "ci"
        ) %>%
            mutate(
                cellularBarcode = gsub(
                    pattern = "-",
                    replacement = "_",
                    x = !!sym("cellularBarcode")
                )
            )
    })
    names(list) <- names(sampleDirs)
    list
}
