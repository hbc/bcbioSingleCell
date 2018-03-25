#' Sanitize Markers Output
#'
#' @note [Seurat::FindAllMarkers()] maps the counts matrix rownames correctly in
#'   the `gene` column, whereas [Seurat::FindMarkers()] maps them correctly in
#'   the rownames of the returned marker `data.frame`.
#'
#' @name sanitizeMarkers
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param markers Original [Seurat::FindAllMarkers()] `data.frame`.
#'
#' @return Tibble arranged by adjusted P value.
#' - [Seurat::FindMarkers()]: `tbl_df`.
#' - [Seurat::FindAllMarkers()]: `grouped_df`, grouped by cluster ident.
#'
#' @examples
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # seurat ====
#' all_markers <- FindAllMarkers(seurat_small)
#' all_sanitized <- sanitizeMarkers(
#'     object = seurat_small,
#'     markers = ident_3_markers
#' )
#' glimpse(all_sanitized)

#' ident_3_markers <- FindMarkers(seurat_small, ident.1 = "3", ident.2 = NULL)
NULL



# Constructors =================================================================
.sanitizeMarkers.seurat <- function(  # nolint
    object,
    markers
) {
    validObject(object)

    # Early return on sanitized data
    if (.isSanitizedMarkers(markers, package = "Seurat")) {
        inform("Markers are already sanitized")
        return(markers)
    }

    assert_has_rows(markers)
    assertHasRownames(markers)

    data <- markers

    # Map the matrix rownames to `rownames` column in tibble
    if ("cluster" %in% colnames(data)) {
        # `FindAllMarkers()` return
        all <- TRUE
        assert_is_subset("gene", colnames(data))
        data[["rowname"]] <- data[["gene"]]
        data[["gene"]] <- NULL
    } else {
        # `FindMarkers()` return
        all <- FALSE
        data <- rownames_to_column(data)
    }

    # Now ready to coerce
    data <- data %>%
        as_tibble() %>%
        # Sanitize column names into lowerCamelCase
        camel()

    # Update legacy columns
    if ("avgDiff" %in% colnames(data)) {
        inform(paste(
            "Renaming legacy `avgDiff` column to `avgLogFC`",
            "(changed in Seurat v2.1)"
        ))
        data[["avgLogFC"]] <- data[["avgDiff"]]
        data[["avgDiff"]] <- NULL
    }

    # Rename P value columns to match DESeq2 conventions
    if ("pVal" %in% colnames(data)) {
        inform("Renaming `pVal` column to `pvalue` (matching DESeq2)")
        data[["pvalue"]] <- data[["pVal"]]
        data[["pVal"]] <- NULL
    }
    if ("pValAdj" %in% colnames(data)) {
        inform("Renaming `pValAdj` column to `padj` (matching DESeq2)")
        data[["padj"]] <- data[["pValAdj"]]
        data[["pValAdj"]] <- NULL
    }

    # Add Ensembl gene IDs
    rownames <- rownames(slot(object, "data"))
    assert_has_names(rownames)
    map <- tibble(
        "rowname" = rownames,
        "geneID" = names(rownames)
    )
    data <- left_join(data, map, by = "rowname")

    # Add genomic ranges, if available
    rowData <- rowData(object)
    if (is.data.frame(rowData)) {
        assert_is_subset("geneID", colnames(rowData))
        # Ensure any nested list columns are dropped
        cols <- vapply(
            X = rowData,
            FUN = function(x) {
                !is.list(x)
            },
            FUN.VALUE = logical(1L)
        )
        rowData <- rowData[, cols, drop = FALSE]
        data <- left_join(data, rowData, by = "geneID")
    }

    # Ensure that required columns are present
    requiredCols <- c(
        "geneID",       # Ensembl annotations
        "geneName",     # Renamed from `gene`
        "pct1",
        "pct2",
        "avgLogFC",     # Seurat v2.1
        "padj",
        "pvalue"        # Renamed from `p_val`
    )
    assert_is_subset(requiredCols, colnames(data))

    if (isTRUE(all)) {
        # `cluster` is only present in `FindAllMarkers() return`
        data <- group_by(data, !!sym("cluster"))
        priorityCols <- c("cluster", "rowname", requiredCols)
        arrangeCols <- c("cluster", "padj")
    } else {
        priorityCols <- c("rowname", requiredCols)
        arrangeCols <- "padj"
    }

    data %>%
        .[, unique(c(priorityCols, colnames(.))), drop = FALSE] %>%
        arrange(!!!syms(arrangeCols))
}



# Methods ======================================================================
#' @rdname sanitizeMarkers
#' @export
setMethod(
    "sanitizeMarkers",
    signature(
        object = "seurat",
        markers = "data.frame"
    ),
    .sanitizeMarkers.seurat
)
