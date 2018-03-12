#' Cell Types per Cluster
#'
#' @note This function only returns the positive markers per cluster.
#'
#' @rdname cellTypesPerCluster
#' @name cellTypesPerCluster
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param min Minimum number of marker genes per cluster.
#' @param max Maximum number of marker genes per cluster.
#'
#' @return [tibble] grouped by cluster, containing the count (`n`) of
#'   significant known makers per cell type.
#'
#' @examples
#' load(system.file("extdata/knownMarkersDetected.rda", package = "bcbioSingleCell"))
#'
#' cellTypesPerCluster(knownMarkersDetected) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr everything group_by select ungroup
#' @importFrom rlang !!! quos
.cellTypesPerCluster <- function(
    object,
    min = 1L,
    max = Inf) {
    if (attr(object, "vars") != "cell") {
        abort("Markers tibble should be grouped by cell")
    }
    requiredCols <- c(
        "avgLogFC",  # Seurat v2.1
        "cell",      # bcbio
        "cluster",   # Seurat
        "geneID",    # bcbio
        "geneName",  # bcbio
        "padj"       # Seurat v2.1
    )
    assert_is_subset(requiredCols, colnames(object))
    groupCols <- syms(c("cluster", "cell"))

    tbl <- object %>%
        ungroup() %>%
        # Use only positive markers for this approach
        filter(.data[["avgLogFC"]] > 0L) %>%
        select(!!!groupCols, everything()) %>%
        group_by(!!!groupCols) %>%
        arrange(.data[["padj"]], .by_group = TRUE) %>%
        summarize(
            n = n(),
            # Genes are arranged by P value
            geneID = toString(.data[["geneID"]]),
            geneName = toString(.data[["geneName"]])
        ) %>%
        group_by(!!sym("cluster")) %>%
        arrange(dplyr::desc(.data[["n"]]), .by_group = TRUE)

    # Apply minimum and maximum gene cutoffs
    if (is.numeric(min) & min > 1L) {
        tbl <- filter(tbl, .data[["n"]] >= min)
    }
    if (is.numeric(max) & max > 1L) {
        tbl <- filter(tbl, .data[["n"]] <= max)

    }
    if (!nrow(tbl)) {
        return(NULL)
    }

    tbl
}



# Methods =====
#' @rdname cellTypesPerCluster
#' @export
setMethod(
    "cellTypesPerCluster",
    signature("grouped_df"),
    .cellTypesPerCluster)
