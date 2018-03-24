#' Sanitize Markers Output
#'
#' @name sanitizeMarkers
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param markers Original [Seurat::FindAllMarkers()] `data.frame`.
#'
#' @return `grouped_df`, grouped by cluster ident.
#'
#' @examples
#' # seurat ====
#' ident_3_markers <- FindMarkers(seurat_small, ident.1 = "3", ident.2 = NULL)
#' sanitizeMarkers(
#'     object = seurat_small,
#'     markers = ident_3_markers
#' ) %>%
#'     glimpse()
NULL



# Constructors =================================================================
.sanitizeMarkers.seurat <- function(  # nolint
    object,
    markers
) {
    if (.checkSanitizedMarkers(markers, package = "Seurat")) {
        inform("Markers are already sanitized")
        return(markers)
    }

    inform("Coercing to tibble")
    markers <- as_tibble(markers)
    markers <- remove_rownames(markers)

    # Sanitize column names into lowerCamelCase
    inform("Converting columns into camelCase")
    markers <- camel(markers)

    # Rename `gene` column to `rowname`
    if ("gene" %in% colnames(markers)) {
        markers[["rowname"]] <- markers[["gene"]]
        markers[["gene"]] <- NULL
    }

    # Update legacy columns
    if ("avgDiff" %in% colnames(markers)) {
        inform(paste(
            "Renaming legacy `avgDiff` column to `avgLogFC`",
            "(changed in Seurat v2.1)"
        ))
        markers[["avgLogFC"]] <- markers[["avgDiff"]]
        markers[["avgDiff"]] <- NULL
    }

    # Rename P value columns to match DESeq2 conventions
    if ("pVal" %in% colnames(markers)) {
        inform("Renaming `pVal` column to `pvalue` (matching DESeq2)")
        markers[["pvalue"]] <- markers[["pVal"]]
        markers[["pVal"]] <- NULL
    }
    if ("pValAdj" %in% colnames(markers)) {
        inform("Renaming `pValAdj` column to `padj` (matching DESeq2)")
        markers[["padj"]] <- markers[["pValAdj"]]
        markers[["pValAdj"]] <- NULL
    }

    # Add Ensembl gene IDs
    rownames <- rownames(slot(object, "data"))
    assert_has_names(rownames)
    symbol2gene <- tibble(
        "rowname" = rownames,
        "geneID" = names(rownames)
    )
    markers <- left_join(markers, symbol2gene, by = "rowname")

    # Add genomic ranges, if available
    rowData <- rowData(object)
    if (is.data.frame(rowData)) {
        inform("Joining row data")
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
        markers <- left_join(markers, rowData, by = "geneID")
    }

    # Ensure that required columns are present
    requiredCols <- c(
        "cluster",      # Unmodified
        "geneID",       # Ensembl annotations
        "geneName",     # Renamed from `gene`
        "pct1",
        "pct2",
        "avgLogFC",     # Seurat v2.1
        "padj",
        "pvalue"        # Renamed from `p_val`
    )
    assert_is_subset(requiredCols, colnames(markers))

    # Return as a tibble
    # Grouped by cluster
    # Arranged by P value (per cluster)
    markers %>%
        camel() %>%
        .[, unique(c(requiredCols, colnames(.))), drop = FALSE] %>%
        group_by(.data[["cluster"]]) %>%
        # Arrange by adjusted P value
        arrange(!!sym("padj"), .by_group = TRUE)
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
