#' Sanitize Markers Output
#'
#' @rdname sanitizeMarkers
#' @name sanitizeMarkers
#'
#' @inheritParams general
#'
#' @param markers [Seurat::FindAllMarkers()] return.
#'
#' @return `grouped_df`, grouped by cluster ident.
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#' load(system.file(
#'     "extdata/seuratAllMarkersOriginal.rda",
#'     package = "bcbioSingleCell"
#' ))
#'
#' # seurat ====
#' sanitizeMarkers(
#'     seurat,
#'     markers = seuratAllMarkersOriginal
#' ) %>%
#'     glimpse()
NULL



# Constructors =================================================================
#' @importFrom basejump camel
#' @importFrom dplyr arrange everything group_by left_join mutate select
#' @importFrom rlang !! sym
#' @importFrom tibble as_tibble remove_rownames
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
    symbol2gene <- tibble(
        "rowname" = rownames,
        "geneID" = names(rownames)
    )
    markers <- left_join(markers, symbol2gene, by = "rowname")

    # Add genomic ranges
    rowRanges <- bcbio(seurat, "rowRanges")
    assert_is_all_of(rowRanges, "GRanges")
    rowData <- as.data.frame(rowRanges)
    markers <- left_join(markers, rowData, by = "geneID")

    # Ensure that required columns are present
    requiredCols <- c(
        "avgLogFC",     # Seurat v2.1
        "cluster",      # Unmodified
        "description",  # Ensembl annotations
        "geneBiotype",  # Ensembl annotations
        "geneName",     # Renamed from `gene`
        "geneID",       # Ensembl annotations
        "pvalue",       # Renamed from `p_val`
        "padj"
    )
    assert_is_subset(requiredCols, colnames(markers))

    # Return as a tibble
    # Grouped by cluster
    # Arranged by P value (per cluster)
    markers %>%
        camel() %>%
        # `padj` should come at the end, but isn't in legacy output
        select(
            c(
                "cluster",
                "geneID",
                "geneName",
                "pct1",
                "pct2",
                "avgLogFC",
                "pvalue",
                "padj"
            ),
            everything()
        ) %>%
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
