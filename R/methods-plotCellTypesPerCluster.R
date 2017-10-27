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
#'
#' @return [ggplot].
#'
#' @examples
#' \dontrun{
#' cellTypes <- cellTypesPerCluster(knownMakersDetected)
#' plotCellTypesPerCluster(cellTypes)
#' }
NULL



# Constructors ====
#' @importFrom dplyr pull
#' @importFrom pbapply pblapply
.plotCellTypesPerCluster <- function(
    object,
    cellTypesPerCluster,
    color = scale_color_viridis(option = "inferno"),
    dark = TRUE,
    headerLevel = NULL) {
    # Output Markdown headers per cluster
    clusters <- unique(cellTypesPerCluster[["cluster"]])
    pblapply(seq_along(clusters), function(a) {
        cluster <- clusters[[a]]
        if (!is.null(headerLevel)) {
            mdHeader(
                paste("Cluster", cluster),
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE)
        }
        subset <- cellTypesPerCluster %>%
                .[.[["cluster"]] == cluster, ]
        lapply(seq_len(nrow(subset)), function(b) {
            cellType <- subset[b, ]
            genes <- pull(cellType, "symbol") %>%
                strsplit(", ") %>%
                .[[1]]
            title <- pull(cellType, "cell")
            if (!is.null(headerLevel)) {
                if (!is.null(headerLevel)) {
                    mdHeader(
                        title,
                        level = headerLevel + 1,
                        tabset = TRUE,
                        asis = TRUE)
                }
            }
            # Modify the title by adding the cluster number (for the plot)
            title <- paste(paste0("Cluster ", cluster, ":"), title)
            plotMarkerTSNE(
                object = object,
                genes = genes,
                colorPoints = "geomean",
                color = color,
                dark = dark,
                pointsAsNumbers = FALSE,
                label = TRUE,
                title = title) %>%
                show()
        }) %>%
            invisible()
    }) %>%
        invisible()
}



# Methods ====
#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    "plotCellTypesPerCluster",
    signature(object = "seurat",
              cellTypesPerCluster = "grouped_df"),
    .plotCellTypesPerCluster)
