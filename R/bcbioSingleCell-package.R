#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @import basejump methods SummarizedExperiment
"_PACKAGE"

globalVariables(".")

bins <- 100L

# Quality control plot colors
qcCutoffColor <- "black"
qcFailColor <- "red"
qcPassColor <- "green"
qcWarnColor <- "orange"

# Quality control label appearance
qcLabelAlpha <- 0.75
qcLabelColor <- "white"
qcLabelFill <- "black"
qcLabelFontface <- "bold"
qcLabelPadding <- ggplot2::unit(0.2, "lines")
qcLabelSize <- NA

# Quality control line appearance
qcLineAlpha <- 0.75
qcLineSize <- 1L
qcLineType <- "dashed"

# Plot label separator
labelSep <- ": "

metaPriorityCols <- c("sampleID", "sampleName", "fileName")
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"

url <- "http://bioinformatics.sph.harvard.edu/bcbioSingleCell"
