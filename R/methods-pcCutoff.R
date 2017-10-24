#' Determine PC Cutoff
#'
#' Calculate the principal component (PC) with either a maximum SD percentage or
#' minimum % SD cumulative sum cutoff value. This defaults to keeping the larger
#' PC cutoff for either 5% SD (`maxPct`) or 80% cumulative (`minCumsum`).
#'
#' @rdname pcCutoff
#' @name pcCutoff
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param maxPct Maximum percent standard deviation.
#' @param minCumPct Minimum cumulative percent standard deviation.
#' @param plot Plot the PC standard deviations.
#'
#' @return Numeric of maximum PC value to use for dimensionality reduction.
#'
#' @seealso [Seurat::PCElbowPlot].
NULL



# Constructors ====
#' @importFrom cowplot draw_plot ggdraw
#' @importFrom ggplot2 aes_string geom_hline geom_line geom_point ggplot labs
#'   scale_y_continuous
#' @importFrom scales percent
#' @importFrom tibble tibble
.pcCutoff <- function(sd, maxPct, minCumPct, plot) {
    xlab <- "pc"

    # Principal component standard deviations
    pct <- sd ^ 2 / sum(sd ^ 2)
    cumsum <- cumsum(pct)
    tbl <- tibble(
        pc = seq_along(sd),
        sdev = sd,
        pct = pct,
        cumsum = cumsum)

    cutoffPct <- tbl %>%
        .[.[["pct"]] >= maxPct, "pc"] %>%
        max()
    cutoffCumPct <- tbl %>%
        .[.[["cumsum"]] <= minCumPct, "pc"] %>%
        max()

    # Pick the larger of the two cutoffs
    cutoff <- max(cutoffPct, cutoffCumPct)

    # Elbow plot
    ggelbow <- ggplot(
        tbl,
        mapping = aes_string(
            x = "pc",
            y = "sd")
    ) +
        geom_point() +
        geom_line() +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = xlab,
             y = "std dev")

    # Percentage plot
    ggpct <- ggplot(
        tbl,
        mapping = aes_string(
            x = "pc",
            y = "pct")
    ) +
        geom_point() +
        geom_line() +
        geom_hline(
            alpha = 0.5,
            color = "orange",
            size = 1.5,
            yintercept = maxPct) +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = xlab,
             y = "% std dev") +
        scale_y_continuous(labels = scales::percent)

    # Cumulative plot
    ggcumsum <- ggplot(
        tbl,
        mapping = aes_string(
            x = "pc",
            y = "cumsum")
    ) +
        geom_point() +
        geom_line() +
        geom_hline(
            alpha = 0.5,
            color = "orange",
            size = 1.5,
            yintercept = minCumPct) +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = xlab,
             y = "cumulative % std dev") +
        scale_y_continuous(labels = scales::percent)

    p <- ggdraw() +
        # Coordinates are relative to lower left corner
        draw_plot(
            ggelbow,
            x = 0, y = 0.5, width = 1, height = 0.5) +
        draw_plot(
            ggpct,
            x = 0, y = 0, width = 0.5, height = 0.5) +
        draw_plot(
            ggcumsum, x = 0.5, y = 0, width = 0.5, height = 0.5)
    show(p)

    cutoff
}



# Methods ====
#' @rdname pcCutoff
#' @export
setMethod("pcCutoff", "seurat", function(
    object,
    maxPct = 0.05,
    minCumPct = 0.8,
    plot = TRUE) {
    # seurat slot descriptions
    # dr: dimensionality reduction
    # sdev: standard deviation
    sd <- slot(object, "dr") %>%
        .[["pca"]] %>%
        slot("sdev")
    .pcCutoff(
        sd,
        maxPct = maxPct,
        minCumPct = minCumPct,
        plot = plot)
})
