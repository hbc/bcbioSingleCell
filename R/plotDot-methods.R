#' Plot Dot
#'
#' @name plotDot
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotDot
#'
#' @inheritParams general
#' @param colMin `scalar numeric`. Minimum scaled average expression threshold.
#'   Everything smaller will be set to this.
#' @param colMax `scalar numeric`.Maximum scaled average expression threshold.
#'   Everything larger will be set to this.
#' @param dotMin `scalar numeric`.The fraction of cells at which to draw the
#'   smallest dot. All cell groups with less than this expressing the given gene
#'   will have no dot drawn.
#' @param dotScale `scalar numeric`.Scale the size of the points, similar to
#'   `cex`.
#'
#' @seealso Modified version of [Seurat::DotPlot()].
#'
#' @return `ggplot`.
#'
#' @examples
#' # SingleCellExperiment ====
#' object <- indrops_small
#' genes <- head(rownames(object))
#' glimpse(genes)
#' plotDot(object, genes = genes)
NULL



# Constructors =================================================================
#' Min Max
#' @seealso [Seurat:::MinMax()].
#' @noRd
.minMax <- function(data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    data2
}



#' Percent Above
#' @seealso [Seurat:::PercentAbove()].
#' @noRd
.percentAbove <- function(x, threshold) {
    length(x[x > threshold]) / length(x)
}



# Methods ======================================================================
#' @rdname plotDot
#' @export
setMethod(
    "plotDot",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        colMin = -2.5,
        colMax = 2.5,
        dotMin = 0L,
        dotScale = 6L,
        color = getOption("bcbio.discrete.color", NULL),
        legend = getOption("bcbio.legend", TRUE)
    ) {
        .assertHasIdent(object)
        assert_is_character(genes)
        assert_is_a_number(colMin)
        assert_is_a_number(colMax)
        assert_is_a_number(dotMin)
        assert_is_a_number(dotScale)
        assertIsColorScaleContinuousOrNULL(color)
        assert_is_a_bool(legend)

        ident <- colData(object)[["ident"]]
        assert_is_non_empty(ident)

        data <- fetchGeneData(
            object = object,
            genes = genes,
            gene2symbol = TRUE
        ) %>%
            as.data.frame() %>%
            cbind(ident) %>%
            rownames_to_column("cell") %>%
            as_tibble()

        if (!isTRUE(.useGene2symbol(object))) {
            g2s <- gene2symbol(object)
            if (length(g2s)) {
                g2s <- g2s[genes, , drop = FALSE]
                genes <- make.unique(g2s[["geneName"]])
                stopifnot(all(genes %in% colnames(data)))
            }
        }

        data <- data %>%
            gather(
                key = "gene",
                value = "expression",
                !!genes
            ) %>%
            group_by(!!!syms(c("ident", "gene"))) %>%
            summarize(
                avgExp = mean(expm1(!!sym("expression"))),
                pctExp = .percentAbove(!!sym("expression"), threshold = 0L)
            ) %>%
            ungroup() %>%
            mutate(gene = factor(!!sym("gene"), levels = genes)) %>%
            group_by(!!sym("gene")) %>%
            mutate(
                avgExpScale = scale(!!sym("avgExp")),
                avgExpScale = .minMax(
                    !!sym("avgExpScale"),
                    max = colMax,
                    min = colMin
                )
            )
        data[["pctExp"]][data[["pctExp"]] < dotMin] <- NA

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("gene"),
                y = !!sym("ident")
            )
        ) +
            geom_point(
                mapping = aes(
                    color = !!sym("avgExpScale"),
                    size = !!sym("pctExp")
                ),
                show.legend = legend
            ) +
            scale_radius(range = c(0L, dotScale)) +
            labs(x = NULL, y = NULL)

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        } else {
            p <- p + lightMarkerColors
        }

        p
    }
)
