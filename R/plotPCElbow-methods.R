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
#' @param minSD `scalar numeric`. Minimum standard deviation.
#' @param minPct `scalar numeric` (`0`-`1`). Minimum percent standard deviation.
#' @param maxCumPct `scalar numeric` (`0`-`1`).Maximum cumulative percent
#'   standard deviation.
#'
#' @return
#' - Show graphical output of elbow plots.
#' - Invisibly return numeric sequence vector of PCs to include for
#'   dimensionality reduction analysis.
#'
#' @seealso
#' - [Seurat::PCElbowPlot()].
#'
#' @examples
#' # seurat ====
#' plotPCElbow(seurat_small)
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
        minPct = 0.01,
        maxCumPct = 0.9,
        trans = c("identity", "sqrt")
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

        # dr: dimensional reduction
        # sdev: standard deviation
        sdev <- object@dr[["pca"]]@sdev
        assert_is_numeric(sdev)
        pct <- sdev ^ 2L / sum(sdev ^ 2L)
        cumsum <- cumsum(pct)

        data <- tibble(
            pc = seq_along(sdev),
            sdev = sdev,
            pct = pct,
            cumsum = cumsum
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
            mapping = aes(
                x = !!sym("pc"),
                y = !!sym("sdev")
            )
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = minSD
            ) +
            geom_line() +
            geom_point() +
            bcbio_geom_abline(xintercept = cutoff) +
            labs(
                x = "pc",
                y = "std dev"
            ) +
            expand_limits(y = 0L) +
            scale_y_continuous(trans = trans)

        # Percent standard deviation ===========================================
        ggpct <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("pc"),
                y = !!sym("pct")
            )
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = minPct
            ) +
            geom_line() +
            geom_point() +
            bcbio_geom_abline(xintercept = cutoff) +
            labs(
                x = "pc",
                y = "% std dev"
            ) +
            expand_limits(y = 0L) +
            scale_y_continuous(labels = percent, trans = trans)

        # Cumulative percent standard deviation ================================
        ggcumsum <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("pc"),
                y = !!sym("cumsum")
            )
        ) +
            geom_hline(
                color = "orange",
                size = 1L,
                yintercept = maxCumPct
            ) +
            geom_line() +
            geom_point() +
            bcbio_geom_abline(xintercept = cutoff) +
            labs(
                x = "pc",
                y = "cum % std dev"
            ) +
            expand_limits(y = c(0L, 1L)) +
            scale_y_continuous(labels = percent, trans = trans)

        plotlist <- list(
            sd = ggsd,
            pct = ggpct,
            cumsum = ggcumsum
        )

        p <- plot_grid(plotlist = plotlist)
        show(p)

        invisible(seq_len(cutoff))
    }
)
