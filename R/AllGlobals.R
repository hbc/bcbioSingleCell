globalVariables(".")

.version <- packageVersion("bcbioSingleCell")

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioSingleCellTestsURL
bcbioSingleCellTestsURL <- paste0(
    "https://tests.acidgenomics.com/bcbioSingleCell/",
    "v", .version$major, ".", .version$minor  # nolint
)

#' @importFrom basejump lanePattern
lanePattern <- basejump::lanePattern

#' @importFrom bcbioBase projectDirPattern
projectDirPattern <- bcbioBase::projectDirPattern

separatorBar <- basejump::separator()

requiredAssays <- "counts"

## This is also defined in acidplots.
geom <- c("histogram", "ecdf", "violin", "ridgeline", "boxplot")

## We're adding an additional raw reads column (pre-UMI disambiguation).
metricsCols <- basejump::metricsCols
metricsCols <- c("nRead", metricsCols)
