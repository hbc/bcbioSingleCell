#' Plot PC Elbow
#'
#' Calculate the principal component (PC) cutoff using a heuristic approach.
#'
#' Automatically return the smallest number of PCs that match the `minSD`,
#' `minPct`, and `maxCumPct` cutoffs.
#'
#' @name plotPCElbow
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param minSD Minimum standard deviation.
#' @param minPct Minimum percent standard deviation.
#' @param maxCumPct Maximum cumulative percent standard deviation.
#' @param plot Include plot.
#' @param trans Name of the transformation to apply. Supports "`identity`"
#'   (unmodified) or "`sqrt`". See [ggplot2::scale_y_continuous()] for more
#'   information.
#'
#' @return
#' - Show graphical output of elbow plots.
#' - Invisibly return numeric sequence vector of PCs to include for
#'   dimensionality reduction analysis.
#'
#' @seealso
#' - [Seurat::PCElbowPlot()].
#' - [ggplot2::scale_y_continuous()].
#'
#' @examples
#' # seurat ====
#' plotPCElbow(pbmc_small)
NULL



# Constructors =================================================================
.plotPCElbow.seurat <- function(  # nolint
    object,
    minSD = 1L,
    minPct = 0.025,
    maxCumPct = 0.9,
    trans = c("identity", "sqrt"),
    plot = TRUE
) {
    assert_is_a_number(minSD)
    assert_all_are_positive(minSD)
    assert_is_a_number(minPct)
    assert_is_a_number(maxCumPct)
    assert_all_are_in_left_open_range(
        x = c(minPct, maxCumPct),
        lower = 0L,
        upper = 1L
    )
    trans <- match.arg(trans)
    assert_is_a_bool(plot)

    # dr: dimensionality reduction
    # sdev: standard deviation
    sdev <- object@dr[["pca"]]@sdev
    assert_is_numeric(sdev)
    pct <- sdev ^ 2L / sum(sdev ^ 2L)
    cumsum <- cumsum(pct)

    data <- tibble(
        "pc" = seq_along(sdev),
        "sdev" = sdev,
        "pct" = pct,
        "cumsum" = cumsum
    )

    minSDCutoff <- data %>%
        .[.[["sdev"]] >= minSD, "pc"] %>%
        max()
    minPctCutoff <- data %>%
        .[.[["pct"]] >= minPct, "pc"] %>%
        max()
    maxCumPctCutoff <- data %>%
        .[.[["cumsum"]] <= maxCumPct, "pc"] %>%
        max()

    # Pick the smallest value of the cutoffs
    cutoff <- min(minSDCutoff, minPctCutoff, maxCumPctCutoff)

    # Standard deviation =======================================================
    ggelbow <- ggplot(
        data = data,
        mapping = aes_string(x = "pc", y = "sdev")
    ) +
        geom_point() +
        geom_line() +
        geom_hline(
            alpha = 0.5,
            color = "orange",
            size = 1.5,
            yintercept = minSD
        ) +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = "pc", y = "std dev") +
        expand_limits(y = 0L) +
        scale_y_continuous(trans = trans)

    # Percent standard deviation ===============================================
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
            yintercept = minPct
        ) +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = "pc", y = "% std dev") +
        expand_limits(y = 0L) +
        scale_y_continuous(labels = percent, trans = trans)

    # Cumulative percent standard deviation ====================================
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
            yintercept = maxCumPct
        ) +
        .qcCutoffLine(xintercept = cutoff) +
        labs(x = "pc", y = "cum % std dev") +
        expand_limits(y = c(0L, 1L)) +
        scale_y_continuous(labels = percent, trans = trans)

    # Coordinates are relative to lower left corner
    p <- ggdraw() +
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
    .plotPCElbow.seurat
)
