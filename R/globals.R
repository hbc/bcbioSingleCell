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
    ## Raw reads (pre UMI disambiguation).
    ## Previously: "nCount" until v0.3.19. Note the switch!
    "nRead",
    ## UMI disambiguated counts.
    ## Previously: "nUMI" until v0.3.19. Note the switch!
    "nCount",
    ## Features (i.e. genes, transcripts).
    ## Previously: "nGene" until v0.3.19.
    ## Rename improves consistency with Chromium code.
    "nFeature",
    ## Coding features.
    "nCoding",
    ## Mitochondrial features.
    "nMito",
    ## Novelty score.
    ## Previously: log10GenesPerUMI until v0.3.19.
    ## Rename improves consistency with Chromium code.
    "log10FeaturesPerCount",
    ## Proportion of mito expression.
    ## Useful to see if we have a lot of stress or dying cells.
    "mitoRatio"
)

Rle <- structure("Rle", package = "S4Vectors")  # nolint
