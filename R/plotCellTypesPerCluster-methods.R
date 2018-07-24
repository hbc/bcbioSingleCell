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
#' @param cellTypesPerCluster `grouped_df`. Cell types per cluster data. This
#'   must be the return from [cellTypesPerCluster()].
#' @param ... Passthrough arguments to [plotMarkerTSNE()] or [plotMarkerUMAP()].
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' # seurat ====
#' per_cluster <- cellTypesPerCluster(known_markers_small)
#' glimpse(per_cluster)
#'
#' # Let's plot the first row, as an example
#' plotCellTypesPerCluster(
#'     object = seurat_small,
#'     cellTypesPerCluster = head(per_cluster, 1),
#' )
NULL



# Methods ======================================================================
#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    "plotCellTypesPerCluster",
    signature("seurat"),
    function(
        object,
        cellTypesPerCluster,
        reduction = c("TSNE", "UMAP"),
        headerLevel = 2L,
        ...
    ) {
        # Passthrough: color, dark
        validObject(object)
        stopifnot(is(cellTypesPerCluster, "grouped_df"))
        assert_has_rows(cellTypesPerCluster)
        assert_are_identical(
            x = group_vars(cellTypesPerCluster),
            y = "cluster"
        )
        reduction <- match.arg(reduction)
        assertIsAHeaderLevel(headerLevel)

        cellTypesPerCluster <- cellTypesPerCluster %>%
            ungroup() %>%
            mutate_if(is.factor, droplevels)

        # Output Markdown headers per cluster
        clusters <- levels(cellTypesPerCluster[["cluster"]])
        assert_is_non_empty(clusters)

        return <- pblapply(clusters, function(cluster) {
            markdownHeader(
                text = paste("Cluster", cluster),
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
                    text = title,
                    level = headerLevel + 1L,
                    asis = TRUE
                )
                # Modify the title by adding the cluster number (for the plot)
                title <- paste(paste0("Cluster ", cluster, ":"), title)
                p <- .plotMarkerReduction(
                    object = object,
                    genes = genes,
                    reduction = reduction,
                    ...
                )
                show(p)
                invisible(p)
            })
        })
        invisible(return)
    }
)
