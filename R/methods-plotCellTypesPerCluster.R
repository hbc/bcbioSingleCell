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
#' @inheritParams AllGenerics
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
#' @importFrom stringr str_split
.plotCellTypesPerCluster <- function(object, cellTypesPerCluster) {
    pblapply(seq_len(nrow(cellTypesPerCluster)), function(a) {
        cellType <- cellTypesPerCluster[a, ]
        genes <- pull(cellType, "symbol") %>%
            str_split(", ") %>%
            .[[1]]
        title <- paste(
            paste("cluster", pull(cellType, "cluster")),
            paste(tolower(pull(cellType, "cell")), "markers"),
            sep = " : "
        )
        plotMarkerTSNE(
            object = object,
            genes = genes,
            colorPoints = "geomean",
            pointsAsNumbers = FALSE,
            label = TRUE,
            title = title) %>%
            show()
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
