#' @importFrom ggplot2 aes_ aes_string coord_flip element_blank element_line
#'   element_rect element_text expand_limits facet_wrap geom_bar geom_boxplot
#'   geom_histogram geom_hline geom_label geom_line geom_point geom_smooth
#'   geom_text geom_violin geom_vline ggplot ggtitle guide_colorbar
#'   guide_legend guides labs qplot scale_color_gradient scale_radius
#'   scale_x_log10 scale_x_sqrt scale_y_continuous scale_y_log10 scale_y_sqrt
#'   theme xlab xlim ylab



#' @importFrom grid unit
#' @importFrom viridis inferno

bins <- 200

lineColor <- "black"
qcColors <- inferno(3)

# Quality control plot colors
qcPassColor <- qcColors[[1]]
qcWarnColor <- qcColors[[2]]
qcFailColor <- qcColors[[3]]
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
        stop("xintercept or yintercept required", call. = FALSE)
    }
}
