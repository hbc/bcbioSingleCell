#' Sanitize Seurat Markers
#'
#' Currently only Seurat markers are supported.
#'
#' @note [Seurat::FindAllMarkers()] maps the counts matrix rownames correctly in
#'   the `gene` column, whereas [Seurat::FindMarkers()] maps them correctly in
#'   the rownames of the returned marker `data.frame`.
#'
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param data `data.frame`. [Seurat::FindAllMarkers()] or
#'   [Seurat::FindMarkers()] return.
#' @param rowRanges `GRanges`. Gene annotations. Must contain `geneName` column
#'   that corresponds to the gene names in the `seurat` object used to generate
#'   the markers `data.frame`.
#'
#' @return `data.frame`, arranged by adjusted P value.
#' @export
#'
#' @examples
#' # seurat ====
#' object <- seurat_small
#'
#' # `FindAllMarkers()` return
#' invisible(capture.output(
#'     all_markers <- Seurat::FindAllMarkers(object)
#' ))
#' all_sanitized <- sanitizeSeuratMarkers(
#'     data = all_markers,
#'     rowRanges = rowRanges(object)
#' )
#' glimpse(all_sanitized)
#'
#' # `FindMarkers()` return
#' invisible(capture.output(
#'     ident_3_markers <- Seurat::FindMarkers(
#'         object = object,
#'         ident.1 = "3",
#'         ident.2 = NULL
#'     )
#' ))
#' ident_3_sanitized <- sanitizeSeuratMarkers(
#'     data = ident_3_markers,
#'     rowRanges = rowRanges(object)
#' )
#' glimpse(ident_3_sanitized)
sanitizeSeuratMarkers <- function(data, rowRanges) {
    assert_is_data.frame(data)
    # Early return on sanitized data
    if (.isSanitizedMarkers(data, package = "Seurat")) {
        message("Markers are already sanitized")
        return(data)
    }
    assertHasRownames(data)
    assert_is_all_of(rowRanges, "GRanges")
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = colnames(mcols(rowRanges))
    )

    # Map the Seurat matrix rownames to `rownames` column in tibble
    if ("cluster" %in% colnames(data)) {
        message("`Seurat::FindAllMarkers()` return detected")
        all <- TRUE
        assert_is_subset("gene", colnames(data))
        data <- remove_rownames(data)
        data[["rowname"]] <- data[["gene"]]
        data[["gene"]] <- NULL
    } else {
        message("`Seurat::FindMarkers()` return detected")
        all <- FALSE
        data <- rownames_to_column(data)
    }

    stopifnot(all(data$rowname %in% rownames(seurat)))
    stopifnot(all(data$rowname %in% names(rowRanges)))

    # Now ready to coerce
    data <- data %>%
        as_tibble() %>%
        # Sanitize column names into lowerCamelCase
        camel()

    # Update legacy columns
    if ("avgDiff" %in% colnames(data)) {
        message(paste(
            "Renaming legacy `avgDiff` column to `avgLogFC`",
            "(changed in Seurat v2.1)"
        ))
        data[["avgLogFC"]] <- data[["avgDiff"]]
        data[["avgDiff"]] <- NULL
    }

    # Rename P value columns to match DESeq2 conventions
    if ("pVal" %in% colnames(data)) {
        data[["pvalue"]] <- data[["pVal"]]
        data[["pVal"]] <- NULL
    }
    if ("pValAdj" %in% colnames(data)) {
        data[["padj"]] <- data[["pValAdj"]]
        data[["pValAdj"]] <- NULL
    }

    # Strip out unwanted seurat columns from rowRanges
    mcols(rowRanges) <- mcols(rowRanges) %>%
        .[!grepl("^gene($|\\.)", colnames(.))]

    # Row data from GRanges
    rowData <- as.data.frame(rowRanges)
    # Ensure any nested list columns are dropped
    cols <- vapply(
        X = rowData,
        FUN = function(x) {
            !is.list(x)
        },
        FUN.VALUE = logical(1L)
    )
    rowData <- rowData[, cols, drop = FALSE]
    rowData[["rowname"]] <- rowData[["geneName"]] %>%
        as.character() %>%
        make.unique()

    # Now safe to join the rowData
    data <- left_join(
        x = data,
        y = rowData,
        by = "rowname"
    )

    # Check that all rows match a geneID
    stopifnot(!any(is.na(data[["geneID"]])))

    # Ensure that required columns are present
    requiredCols <- c(
        "rowname",
        "geneID",
        "geneName",
        "pct1",
        "pct2",
        "avgLogFC",     # Seurat v2.1
        "padj",
        "pvalue"        # Renamed from `p_val`
    )
    assert_is_subset(requiredCols, colnames(data))

    if (isTRUE(all)) {
        # `cluster` is only present in `FindAllMarkers() return`
        data <- data %>%
            select(!!sym("cluster"), everything()) %>%
            group_by(!!sym("cluster")) %>%
            arrange(!!sym("padj"), .by_group = TRUE)
    } else {
        data <- data %>%
            arrange(!!sym("padj")) %>%
            as.data.frame() %>%
            column_to_rownames()
    }

    message("Sanitized to contain `geneID` and `geneName` columns from GRanges")
    data
}
