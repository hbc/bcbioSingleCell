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
#' plotPCElbow(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotPCElbow
#' @export
setMethod(
    "plotPCElbow",
    signature("seurat"),
    function(
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

        # Standard deviation ===================================================
        ggsd <- ggplot(
            data = data,
            mapping = aes_string(x = "pc", y = "sdev")
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = minSD
            ) +
            geom_line() +
            geom_point() +
            bcbio_geom_abline(xintercept = cutoff) +
            labs(x = "pc", y = "std dev") +
            expand_limits(y = 0L) +
            scale_y_continuous(trans = trans)

        # Percent standard deviation ===========================================
        ggpct <- ggplot(
            data = data,
            mapping = aes_string(x = "pc", y = "pct")
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = minPct
            ) +
            geom_line() +
            geom_point() +
            bcbio_geom_abline(xintercept = cutoff) +
            labs(x = "pc", y = "% std dev") +
            expand_limits(y = 0L) +
            scale_y_continuous(labels = percent, trans = trans)

        # Cumulative percent standard deviation ================================
        ggcumsum <- ggplot(
            data = data,
            mapping = aes_string(x = "pc", y = "cumsum")
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = maxCumPct
            ) +
            geom_line() +
            geom_point() +
            bcbio_geom_abline(xintercept = cutoff) +
            labs(x = "pc", y = "cum % std dev") +
            expand_limits(y = c(0L, 1L)) +
            scale_y_continuous(labels = percent, trans = trans)

        plotlist <- list(
            sd = ggsd,
            pct = ggpct,
            cumsum = ggcumsum
        )

        p <- plot_grid(plotlist = plotlist, labels = "AUTO")
        show(p)

        invisible(seq_len(cutoff))
    }
)
