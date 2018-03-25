#' Cell Types per Cluster
#'
#' @note This function only returns the positive markers per cluster.
#'
#' @name cellTypesPerCluster
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param min Minimum number of marker genes per cluster.
#' @param max Maximum number of marker genes per cluster.
#'
#' @return `tibble` grouped by cluster, containing the count (`n`) of
#'   significant known makers per cell type.
#'
#' @examples
#' load(system.file(
#'     "extdata/known_markers_small.rda",
#'     package = "bcbioSingleCell"
#' ))
#' cellTypesPerCluster(known_markers_small) %>% glimpse()
NULL



# Constructors =================================================================
.cellTypesPerCluster <- function(
    object,
    min = 1L,
    max = Inf
) {
    assert_are_identical(attr(object, "vars"), "cellType")
    requiredCols <- c(
        "cellType",  # bcbio
        "cluster",   # Seurat
        "geneID",    # bcbio
        "geneName",  # bcbio
        "avgLogFC",  # Seurat v2.1
        "padj"       # Seurat v2.1
    )
    assert_is_subset(requiredCols, colnames(object))

    # Note that the order is important here
    groupCols <- c("cluster", "cellType")

    data <- object %>%
        ungroup() %>%
        .[, unique(c(groupCols, colnames(.))), drop = FALSE] %>%
        group_by(!!!syms(groupCols)) %>%
        arrange(.data[["padj"]], .by_group = TRUE) %>%
        # Use `toString()` instead of `aggregate()` for R Markdown tables
        summarize(
            n = n(),
            # Genes are arranged by P value
            geneID = toString(.data[["geneID"]]),
            geneName = toString(.data[["geneName"]])
        ) %>%
        group_by(!!sym("cluster")) %>%
        arrange(dplyr::desc(.data[["n"]]), .by_group = TRUE)

    # Apply minimum and maximum gene cutoffs
    if (is.numeric(min) && min > 1L) {
        data <- data[data[["n"]] >= min, , drop = FALSE]
    }
    if (is.numeric(max) && max > 1L) {
        data <- data[data[["n"]] <= max, , drop = FALSE]
    }
    assert_has_rows(data)

    data
}



# Methods =====
#' @rdname cellTypesPerCluster
#' @export
setMethod(
    "cellTypesPerCluster",
    signature("grouped_df"),
    .cellTypesPerCluster
)
