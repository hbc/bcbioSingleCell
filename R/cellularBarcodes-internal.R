#' Bind Cellular Barcodes
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param list `list`. Cellular barcodes per sample.
#'
#' @return `data.frame`.
.bindCellularBarcodes <- function(list) {
    assert_is_list(list)
    mapply(
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
        select(!!!syms(c("sampleID", "cellularBarcode", "nCount")))
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
    list <- lapply(files, function(file) {
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
