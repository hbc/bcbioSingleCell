# FIXME Simply gate at a standard deviation of 1 cutoff

#' Plot PC Elbow
#'
#' Calculate the principal component (PC) with either a maximum SD percentage or
#' minimum % SD cumulative sum cutoff value. This defaults to keeping the larger
#' PC cutoff for either 5% SD (`maxPct`) or 90% cumulative (`minCumsum`).
#'
#' @name plotPCElbow
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param maxPct Maximum percent standard deviation.
#' @param minCumPct Minimum cumulative percent standard deviation.
#' @param plot Plot the PC standard deviations.
#'
#' @return
#' - Show graphical output of elbow plots.
#' - Invisibly return numeric sequence vector of PCs to include for
#'   dimensionality reduction analysis.
#'
#' @seealso [Seurat::PCElbowPlot()].
#'
#' @examples
#' # seurat ====
#' plotPCElbow(pbmc_small)
NULL



# Constructors =================================================================
.plotPCElbow <- function(
    sd,
    maxPct = 0.05,
    minCumPct = 0.8,
    plot = TRUE
) {
    xlab <- "pc"

    # Principal component standard deviations
    pct <- sd ^ 2L / sum(sd ^ 2L)
    cumsum <- cumsum(pct)
    data <- tibble(
        "pc" = seq_along(sd),
        "sdev" = sd,
        "pct" = pct,
        "cumsum" = cumsum
    )

    cutoffPct <- data %>%
        .[.[["pct"]] >= maxPct, "pc"] %>%
        max()
    cutoffCumPct <- data %>%
        .[.[["cumsum"]] <= minCumPct, "pc"] %>%
        max()

    # Pick the larger of the two cutoffs
    cutoff <- max(cutoffPct, cutoffCumPct)

    # Elbow plot
    ggelbow <- ggplot(
        data = data,
        mapping = aes_string(x = "pc", y = "sd")
    ) +
        geom_point() +
        geom_line() +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = xlab, y = "std dev")

    # Percentage plot
    ggpct <- ggplot(
        data = data,
        mapping = aes_string(x = "pc", y = "pct")
    ) +
        geom_point() +
        geom_line() +
        geom_hline(
            alpha = 0.5,
            color = "orange",
            size = 1.5,
            yintercept = maxPct
        ) +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = xlab, y = "% std dev") +
        scale_y_continuous(labels = percent)

    # Cumulative plot
    ggcumsum <- ggplot(
        data = data,
        mapping = aes_string(x = "pc", y = "cumsum")
    ) +
        geom_point() +
        geom_line() +
        geom_hline(
            alpha = 0.5,
            color = "orange",
            size = 1.5,
            yintercept = minCumPct
        ) +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = xlab, y = "cumulative % std dev") +
        scale_y_continuous(labels = percent)

    p <- ggdraw() +
        # Coordinates are relative to lower left corner
        draw_plot(
            plot = ggelbow,
            x = 0L,
            y = 0.5,
            width = 1L,
            height = 0.5
        ) +
        draw_plot(
            plot = ggpct,
            x = 0L,
            y = 0L,
            width = 0.5,
            height = 0.5
        ) +
        draw_plot(
            plot = ggcumsum,
            x = 0.5,
            y = 0L,
            width = 0.5,
            height = 0.5
        )
    show(p)

    invisible(seq_len(cutoff))
}



# Methods ======================================================================
#' @rdname plotPCElbow
#' @export
setMethod(
    "plotPCElbow",
    signature("seurat"),
    function(
        object,
        maxPct = 0.05,
        minCumPct = 0.9,
        plot = TRUE
    ) {
        # seurat slot descriptions
        # dr: dimensionality reduction
        # sdev: standard deviation
        sd <- object %>%
            slot("dr") %>%
            .[["pca"]] %>%
            slot("sdev")
        .plotPCElbow(
            sd,
            maxPct = maxPct,
            minCumPct = minCumPct,
            plot = plot
        )
    }
)
