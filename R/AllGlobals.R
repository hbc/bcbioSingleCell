## This is also defined in AcidPlots.
.geom <- c("histogram", "ecdf", "violin", "ridgeline", "boxplot")

## We're adding an additional raw reads column (pre-UMI disambiguation).
.metricsCols <- c("nRead", metricsCols)

.requiredAssays <- "counts"

.version <- packageVersion(packageName())

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioSingleCellTestsUrl
bcbioSingleCellTestsUrl <- paste0(
    "https://r.acidgenomics.com/testdata/bcbiosinglecell/",
    "v", .version$major, ".", .version$minor # nolint
)
