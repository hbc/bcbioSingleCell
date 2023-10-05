.pkgName <- packageName()
.pkgVersion <- packageVersion(.pkgName)

## This is also defined in AcidPlots.
.geom <- c("histogram", "ecdf", "violin", "ridgeline", "boxplot")

## We're adding an additional raw reads column (pre-UMI disambiguation).
.metricsCols <- c("nRead", metricsCols)

.requiredAssays <- "counts"

#' Cache URL
#' @keywords internal
#' @export
#' @examples
#' bcbioSingleCellTestsUrl
bcbioSingleCellTestsUrl <- "https://r.acidgenomics.com/testdata/bcbiosinglecell"
