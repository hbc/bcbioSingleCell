bins <- 200L

lineColor <- "black"
qcColors <- viridis::inferno(3L)

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
qcLabelPadding <- grid::unit(0.2, "lines")
qcLabelSize <- NA

# Maximum number of samples to label in a bar or boxplot
qcLabelMaxNum <- 16L

# Quality control line appearance
qcLineAlpha <- 0.75
qcLineSize <- 1L
qcLineType <- "dashed"

# Plot label separator
labelSep <- ": "

inflectionColor <- "red"
kneeColor <- "orange"



# Internal functions ===========================================================
.minimalAxes <- function() {
    theme(
        axis.line = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        title = element_blank()
    )
}



# TODO Migrate this shared code to bcbioBase
.geomLabel <- function(
    data = NULL,
    mapping = NULL,
    color = NULL,
    size = 4L
) {
    geom <- geom_label_repel(
        data = data,
        mapping = mapping,
        arrow = arrow(length = unit(0.01, "npc")),
        box.padding = unit(0.5, "lines"),
        fill = "white",
        fontface = "bold",
        force = 1L,
        point.padding = unit(0.75, "lines"),
        segment.size = 0.5,
        show.legend = FALSE,
        size = size
    )
    if (is.character(color)) {
        geom[["aes_params"]][["colour"]] <- color
    }
    geom
}



# TODO Migrate this shared code to bcbioBase
.qcCutoffLine <- function(
    xintercept = NULL,
    yintercept = NULL,
    alpha = qcLineAlpha,
    color = qcCutoffColor,
    linetype = qcLineType,
    size = qcLineSize
) {
    if (is.null(xintercept) && is.null(yintercept)) {
        stop("`xintercept` and `yintercept` are both NULL")
    } else if (is.numeric(xintercept) && is.numeric(yintercept)) {
        stop("Specifcly only `xintercept` or `yintercept` as numeric")
    } else if (is.numeric(xintercept)) {
        geom_vline(
            xintercept = xintercept,
            alpha = alpha,
            color = color,
            linetype = linetype,
            size = size
        )
    } else if (is.numeric(yintercept)) {
        geom_hline(
            yintercept = yintercept,
            alpha = alpha,
            color = color,
            linetype = linetype,
            size = size
        )
    }
}
