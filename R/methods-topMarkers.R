#' Top Markers
#'
#' @rdname topMarkers
#' @name topMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @param n Number of genes per cluster.
#' @param coding Only include protein coding genes.
#' @param show Show [kable].
#'
#' @seealso [dplyr::top_n()].
#'
#' @return [tibble].
#' @export
NULL



# Constructors ====
# Also called by `knownMarkersDetected()`
.groupMarkers <- function(df) {
    df <- camel(df)

    # [Seurat::FindAllMarkers()] output
    if (identical(colnames(df),
        c("pVal", "avgDiff", "pct1", "pct2", "cluster", "gene"))) {
        # Rename the gene symbol and P value columns
        message("Fixing columns in Seurat marker data.frame")
        df <- df %>%
            rename(pvalue = .data[["pVal"]],
                   symbol = .data[["gene"]])
    }

    # Ensure that required columns are present
    requiredCols <- c("cluster", "ensgene", "symbol", "pvalue")
    if (!all(c("cluster", "ensgene", "symbol", "pvalue") %in% colnames(df))) {
        stop(paste("Marker data.frame must contain:", toString(requiredCols)))
    }

    df %>%
        remove_rownames %>%
        as("tibble") %>%
        camel %>%
        tidy_select(c("cluster", "symbol"), everything()) %>%
        group_by(.data[["cluster"]])
}



# Also called by `knownMarkersDetected()`
.markersKable <- function(object, caption = NULL) {
    object %>%
        .[, c("cluster",
              "symbol",
              "ensgene",
              "pvalue",
              "description",
              "biotype")] %>%
        # Format the P values into consistent scientific notation
        mutate(pvalue = format(.data[["pvalue"]],
                               digits = 3L,
                               scientific = TRUE)) %>%
        # Remove the source information from description
        mutate(description = str_replace(
            .data[["description"]], " \\[.+$", "")) %>%
        # Truncate the description to 50 characters
        mutate(description = str_trunc(.data[["description"]], 50L)) %>%
        # Set digits to a large value to prevent rounding of P values
        # This works but is hacky. I'd like 0.XXXe-XX.
        kable(caption = caption) %>%
        show
}



.topMarkers <- function(
    object,
    n = 4L,
    coding = FALSE,
    show = FALSE) {
    markers <- object
    if (isTRUE(coding)) {
        markers <- tidy_filter(markers, .data[["biotype"]] == "protein_coding")
    }
    markers <- .groupMarkers(markers) %>%
        # Use only the positive markers
        tidy_filter(.data[["avgDiff"]] > 0L) %>%
        # `-n` here means take the smallest P values
        top_n(n = -n, wt = .data[["pvalue"]])
    if (isTRUE(show)) {
        .markersKable(markers, caption = paste("Top", n, "markers per cluster"))
    }
    markers
}



# Methods ====
#' @rdname topMarkers
#' @export
setMethod("topMarkers", "data.frame", .topMarkers)
