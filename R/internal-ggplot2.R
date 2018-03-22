bins <- 200L

lineColor <- "black"
qcColors <- inferno(3L)

# Quality control plot colors
qcPassColor <- qcColors[[1L]]
qcWarnColor <- qcColors[[2L]]
qcFailColor <- qcColors[[3L]]
qcCutoffColor <- qcPassColor

qcPlotAlpha <- 0.85

# Quality control label appearance
qcLabelAlpha <- 0.75
qcLabelColor <- "white"
qcLabelFill <- "black"
qcLabelFontface <- "bold"
qcLabelPadding <- unit(0.2, "lines")
qcLabelSize <- NA

# Maximum number of samples to label in a bar or boxplot
qcLabelMaxNum <- 16L

# Quality control line appearance
qcLineAlpha <- 0.75
qcLineSize <- 1L
qcLineType <- "dashed"

# Plot label separator
labelSep <- ": "



# Internal functions ===========================================================
.qcCutoffLine <- function(xintercept, yintercept) {
    if (!missing(xintercept)) {
        geom_vline(
            alpha = qcLineAlpha,
            color = qcCutoffColor,
            linetype = qcLineType,
            size = qcLineSize,
            xintercept = xintercept
        )
    } else if (!missing(yintercept)) {
        geom_hline(
            alpha = qcLineAlpha,
            color = qcCutoffColor,
            linetype = qcLineType,
            size = qcLineSize,
            yintercept = yintercept
        )
    } else {
        abort("`xintercept` or `yintercept` value required")
    }
}
