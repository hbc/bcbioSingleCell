#' Plot Cell Types per Cluster
#'
#' Plot the geometric mean of the significant marker genes for every known cell
#' type (per unbiased cluster). Cell types with too few (`min` cutoff) or too
#' many (`max` cutoff) marker genes will be skipped.
#'
#' @name plotCellTypesPerCluster
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param cellTypesPerCluster Cell types per cluster `grouped_df`. This must be
#'   the return from [cellTypesPerCluster()].
#'
#' @return Show graphical output. Invisibly return `ggplot` plotlist.
#'
#' @examples
#' # seurat ====
#' per_cluster <- cellTypesPerCluster(known_markers_small)
#' glimpse(per_cluster)
#'
#' # Let's plot the first row, as an example
#' plotCellTypesPerCluster(
#'     object = seurat_small,
#'     cellTypesPerCluster = per_cluster[1L, , drop = FALSE]
#' )
NULL



# Constructors =================================================================
.plotCellTypesPerCluster <- function(
    object,
    cellTypesPerCluster,
    color = scale_color_viridis(discrete = FALSE),
    dark = TRUE,
    headerLevel = 2L
) {
    assert_has_rows(cellTypesPerCluster)
    assert_are_identical(
        x = group_vars(cellTypesPerCluster),
        y = "cluster"
    )
    assertIsColorScaleContinuousOrNULL(color)
    assert_is_a_bool(dark)
    assertIsAHeaderLevel(headerLevel)

    cellTypesPerCluster <- cellTypesPerCluster %>%
        ungroup() %>%
        mutate_if(is.factor, droplevels)

    # Output Markdown headers per cluster
    clusters <- levels(cellTypesPerCluster[["cluster"]])
    assert_is_non_empty(clusters)

    return <- pblapply(clusters, function(cluster) {
        markdownHeader(
            paste("Cluster", cluster),
            level = headerLevel,
            tabset = TRUE,
            asis = TRUE
        )
        subset <- cellTypesPerCluster %>%
            .[.[["cluster"]] == cluster, , drop = FALSE]
        assert_has_rows(subset)
        lapply(seq_len(nrow(subset)), function(x) {
            cellType <- subset[x, , drop = FALSE]
            genes <- pull(cellType, "geneName") %>%
                strsplit(", ") %>%
                .[[1L]]
            title <- pull(cellType, "cellType")
            markdownHeader(
                title,
                level = headerLevel + 1L,
                tabset = TRUE,
                asis = TRUE
            )
            # Modify the title by adding the cluster number (for the plot)
            title <- paste(paste0("Cluster ", cluster, ":"), title)
            p <- plotMarkerTSNE(
                object = object,
                genes = genes,
                expression = "mean",
                color = color,
                dark = dark,
                pointsAsNumbers = FALSE,
                label = TRUE,
                title = title
            )
            show(p)
            invisible(p)
        })
    })
    invisible(return)
}



# Methods ======================================================================
#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    "plotCellTypesPerCluster",
    signature(
        object = "seurat",
        cellTypesPerCluster = "grouped_df"
    ),
    .plotCellTypesPerCluster
)
