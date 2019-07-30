globalVariables(".")

.version <- packageVersion("bcbioSingleCell")

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioSingleCellTestsURL
bcbioSingleCellTestsURL <- paste0(
    "http://tests.acidgenomics.com/bcbioSingleCell/",
    "v", .version$major, ".", .version$minor  # nolint
)

#' @importFrom basejump lanePattern
lanePattern <- basejump::lanePattern

#' @importFrom bcbioBase projectDirPattern
projectDirPattern <- bcbioBase::projectDirPattern

separatorBar <- basejump::separator()

requiredAssays <- "counts"

geom <- c("histogram", "ecdf", "violin", "ridgeline", "boxplot")

## `nCount` column is bcbioSingleCell class specific.
metricsCols <- c(
    "nUMI",
    "nGene",
    "nCoding",
    "nMito",
    "log10GenesPerUMI",
    "mitoRatio"
)

Rle <- structure("Rle", package = "S4Vectors")  # nolint
