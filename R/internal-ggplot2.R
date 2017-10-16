#' @importFrom ggplot2 unit
#' @importFrom viridis inferno

bins <- 200

lineColor <- "black"
qcColors <- viridis::inferno(3)

# Quality control plot colors
qcPassColor <- qcColors[[1]]
qcWarnColor <- qcColors[[2]]
qcFailColor <- qcColors[[3]]
qcCutoffColor <- qcPassColor

qcPlotAlpha <- 0.85
qcRidgeScale <- 2

# Quality control label appearance
qcLabelAlpha <- 0.75
qcLabelColor <- "white"
qcLabelFill <- "black"
qcLabelFontface <- "bold"
qcLabelPadding <- ggplot2::unit(0.2, "lines")
qcLabelSize <- NA

# Maximum number of samples to label in a bar or boxplot
qcLabelMaxNum <- 16

# Quality control line appearance
qcLineAlpha <- 0.75
qcLineSize <- 1
qcLineType <- "dashed"

# Plot label separator
labelSep <- ": "



# Internal functions ====
#' @importFrom ggplot2 geom_hline geom_vline
.qcCutoffLine <- function(xintercept, yintercept) {
    if (!missing(xintercept)) {
        geom_vline(
            alpha = qcLineAlpha,
            color = qcCutoffColor,
            linetype = qcLineType,
            size = qcLineSize,
            xintercept = xintercept)
    } else if (!missing(yintercept)) {
        geom_hline(
            alpha = qcLineAlpha,
            color = qcCutoffColor,
            linetype = qcLineType,
            size = qcLineSize,
            yintercept = yintercept)
    } else {
        stop("xintercept or yintercept required")
    }
}
