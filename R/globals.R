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

## Previously: "nGene", "log10GenesPerUMI" (until v0.3.19).
metricsCols <- c(
    "nCount",
    "nUMI",
    "nFeature",
    "nCoding",
    "nMito",
    "log10FeaturesPerUMI",
    "mitoRatio"
)

Rle <- structure("Rle", package = "S4Vectors")  # nolint
