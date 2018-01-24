#' Plot Cell Types per Cluster
#'
#' Plot the geometric mean of the significant marker genes for every known cell
#' type (per unbiased cluster). Cell types with too few (`min` cutoff) or too
#' many (`max` cutoff) marker genes will be skipped.
#'
#' @rdname plotCellTypesPerCluster
#' @name plotCellTypesPerCluster
#' @author Michael Steinbaugh
#'
#' @inherit plotMarkers
#'
#' @param cellTypesPerCluster Cell types per cluster grouped [tibble]. This must
#'   be the return from [cellTypesPerCluster()].
#' @param color Color palette.
#'
#' @return Show graphical output. Invisibly return [ggplot] plotlist.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "knownMarkersDetected.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' cellTypesPerCluster <- cellTypesPerCluster(knownMarkersDetected)
#' glimpse(cellTypesPerCluster)
#'
#' # seurat
#' # Let's plot the first row, as an example
#' cellTypesPerCluster <- cellTypesPerCluster[1, , drop = FALSE]
#' plotCellTypesPerCluster(
#'     seurat,
#'     cellTypesPerCluster = cellTypesPerCluster)
NULL



# Constructors =================================================================
#' @importFrom dplyr group_vars mutate_if pull
#' @importFrom pbapply pblapply
.plotCellTypesPerCluster <- function(
    object,
    cellTypesPerCluster,
    color = viridis::scale_color_viridis(),
    dark = TRUE,
    headerLevel = NULL) {
    if (!nrow(cellTypesPerCluster)) return(NULL)
    if (group_vars(cellTypesPerCluster) != "cluster") {
        abort("cellTypesPerCluster must be grouped by `cluster` column")
    }
    cellTypesPerCluster <- cellTypesPerCluster %>%
        ungroup() %>%
        mutate_if(is.factor, droplevels)
    # Output Markdown headers per cluster
    clusters <- levels(cellTypesPerCluster[["cluster"]])
    if (is.null(clusters)) return(NULL)
    return <- pblapply(seq_along(clusters), function(a) {
        cluster <- clusters[[a]]
        if (!is.null(headerLevel)) {
            mdHeader(
                paste("Cluster", cluster),
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE)
        }
        subset <- cellTypesPerCluster %>%
                .[.[["cluster"]] == cluster, , drop = FALSE]
        lapply(seq_len(nrow(subset)), function(b) {
            cellType <- subset[b, , drop = FALSE]
            genes <- pull(cellType, "symbol") %>%
                strsplit(", ") %>%
                .[[1L]]
            title <- pull(cellType, "cell")
            if (!is.null(headerLevel)) {
                mdHeader(
                    title,
                    level = headerLevel + 1L,
                    tabset = TRUE,
                    asis = TRUE)
            }
            # Modify the title by adding the cluster number (for the plot)
            title <- paste(paste0("Cluster ", cluster, ":"), title)
            p <- plotMarkerTSNE(
                object = object,
                genes = genes,
                colorPoints = "geomean",
                color = color,
                dark = dark,
                pointsAsNumbers = FALSE,
                label = TRUE,
                title = title)
            show(p)
            p
        })
    })
    invisible(return)
}



# Methods ======================================================================
#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    "plotCellTypesPerCluster",
    signature(object = "seurat",
              cellTypesPerCluster = "grouped_df"),
    .plotCellTypesPerCluster)
