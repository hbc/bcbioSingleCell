#' Plot Dot
#'
#' @name plotDot
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotDot
#'
#' @inheritParams general
#' @param colMin Minimum scaled average expression threshold. Everything
#'   smaller will be set to this.
#' @param colMax Maximum scaled average expression threshold. Everything larger
#'   will be set to this.
#' @param dotMin The fraction of cells at which to draw the smallest dot. All
#'   cell groups with less than this expressing the given gene will have no dot
#'   drawn.
#' @param dotScale Scale the size of the points, similar to `cex`.
#'
#' @seealso Modified version of [Seurat::DotPlot()].
#'
#' @examples
#' # seurat ====
#' genes <- head(rownames(seurat_small))
#' plotDot(seurat_small, genes = genes)
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
    signature("seurat"),
    function(
        object,
        genes,
        color = "auto",
        dark = FALSE,
        grid = FALSE,
        colMin = -2.5,
        colMax = 2.5,
        dotMin = 0L,
        dotScale = 6L
    ) {
        assert_is_character(genes)
        assert_is_a_number(colMin)
        assert_is_a_number(colMax)
        assert_is_a_number(dotMin)
        assert_is_a_number(dotScale)

        ident <- slot(object, "ident")
        data <- fetchGeneData(object, genes = genes) %>%
            as.data.frame() %>%
            cbind(ident) %>%
            rownames_to_column("cell") %>%
            as_tibble() %>%
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
            mapping = aes_string(x = "gene", y = "ident")
        ) +
            geom_point(
                mapping = aes_string(color = "avgExpScale", size = "pctExp"),
                show.legend = FALSE
            ) +
            scale_radius(range = c(0L, dotScale)) +
            labs(x = "gene", y = "ident")

        if (isTRUE(dark)) {
            p <- p + theme_midnight(grid = grid)
            if (color == "auto") {
                color <- scale_color_gradient(
                    low = "gray10",
                    high = "white"
                )
            }
        } else {
            p <- p + theme_paperwhite(grid = grid)
            if (color == "auto") {
                color <- scale_color_gradient(
                    low = "gray90",
                    high = "black"
                )
            }
        }

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        p
    }
)
