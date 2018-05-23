#' Plot Gene
#'
#' @name plotDot
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotGene
#'
#' @param geom Plot type. Uses [match.arg()] to pick the type. Currently
#'   supports "`dot`" and "`violin`".
#'
#' @seealso
#' - [Seurat::DotPlot()].
#' - [Seurat::RidgePlot()].
#' - [Seurat::ViolinPlot()].
#'
#' @examples
#' object <- seurat_small
#' genes <- head(rownames(object))
#'
#' # Dot
#' plotGene(object, genes = genes, geom = "dot")
#'
#' # Violin
#' plotGene(object, genes = genes, geom = "violin")
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



#' Plot Dot
#' @param colMin Minimum scaled average expression threshold. Everything
#'   smaller will be set to this.
#' @param colMax Maximum scaled average expression threshold. Everything larger
#'   will be set to this.
#' @param dotMin The fraction of cells at which to draw the smallest dot. All
#'   cell groups with less than this expressing the given gene will have no dot
#'   drawn.
#' @param dotScale Scale the size of the points, similar to `cex`.
#' @noRd
.plotDot <- function(
    object,
    genes,
    color = "auto",
    dark = FALSE,
    grid = FALSE,
    colMin = -2.5,
    colMax = 2.5,
    dotMin = 0L,
    dotScale = 6L,
    legend = FALSE
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
            show.legend = legend
        ) +
        scale_radius(range = c(0L, dotScale)) +
        labs(x = NULL, y = NULL)

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



.plotViolin <- function(
    object,
    genes,
    scale = c("count", "width", "area"),
    fill = NULL,
    dark = FALSE,
    legend = FALSE
) {
    scale <- match.arg(scale)
    assert_is_any_of(fill, c("ScaleDiscrete", "character", "NULL"))
    if (is.character(fill)) {
        assert_is_a_string(fill)
    }
    assert_is_a_bool(dark)

    if (isTRUE(dark)) {
        color <- "white"
    } else {
        color <- "black"
    }

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
        group_by(!!sym("gene"))

    violin <- geom_violin(
        mapping = aes_string(fill = "ident"),
        # never include a color border
        color = color,
        scale = scale,
        adjust = 1L,
        show.legend = legend,
        trim = TRUE
    )
    if (is_a_string(fill)) {
        violin[["aes_params"]][["fill"]] <- fill
    }

    p <- ggplot(
        data,
        mapping = aes_string(
            x = "ident",
            y = "expression"
        )
    ) +
        violin +
        # Note that `scales = free_y` will hide the x-axis for some plots
        facet_wrap(facets = "gene", scales = "free")

    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    if (isTRUE(dark)) {
        p <- p + theme_midnight()
    } else {
        p <- p + theme_paperwhite()
    }

    p
}



# Methods ======================================================================
#' @rdname plotDot
#' @export
setMethod(
    "plotGene",
    signature("seurat"),
    function(
        object,
        genes,
        geom = c("dot", "violin"),
        legend = FALSE
    ) {
        geom <- match.arg(geom)
        if (geom == "dot") {
            fun <- .plotDot
        } else if (geom == "violin") {
            fun <- .plotViolin
        }
        fun(
            object = object,
            genes = genes,
            legend = legend
        )
    }
)
