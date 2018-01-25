#' Plot Dot
#'
#' @rdname plotDot
#' @name plotDot
#'
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotDot
#'
#' @inheritParams AllGenerics
#'
#' @param genes Gene identifiers to plot.
#' @param format Gene identifier format. Supports `ensgene` or `symbol`.
#' @param color Color palette. If `NULL`, uses default ggplot2 colors.
#' @param colMin Minimum scaled average expression threshold. Everything
#'   smaller will be set to this.
#' @param colMax Maximum scaled average expression threshold. Everything larger
#'   will be set to this.
#' @param dotMin The fraction of cells at which to draw the smallest dot. All
#'   cell groups with less than this expressing the given gene will have no dot
#'   drawn.
#' @param dotScale Scale the size of the points, similar to `cex`.
#'
#' @note [Seurat::DotPlot()] is still plotting even when `do.return = TRUE`.
#' In the meantime, we've broken out the code into this generic to fix
#' RMarkdown looping of our marker plots.
#'
#' @seealso Modified version of [Seurat::DotPlot()].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' symbol <- slot(seurat, "data") %>% rownames() %>% .[1:2]
#' print(symbol)
#'
#' ensgene <- bcbio(seurat, "gene2symbol") %>%
#'     .[which(.[["symbol"]] %in% symbol), "ensgene", drop = TRUE]
#'
#' # seurat
#' plotDot(seurat, genes = symbol, format = "symbol")
#' plotDot(seurat, genes = ensgene, format = "ensgene")
NULL



# Constructors =================================================================
#' Min Max
#' @seealso `Seurat:::MinMax()`
#' @noRd
.minMax <- function(data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    data2
}



#' Percent Above
#' @seealso `Seurat:::PercentAbove()`
#' @noRd
.percentAbove <- function(x, threshold) {
    length(x[x > threshold]) / length(x)
}



#' @importFrom dplyr group_by mutate summarize ungroup
#' @importFrom ggplot2 aes_string geom_point labs scale_color_gradient
#'   scale_radius
#' @importFrom tibble rownames_to_column
#' @importFrom rlang !! !!! sym syms
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
.plotDot.seurat <- function(  # nolint
    object,
    genes,
    format = "symbol",
    color = ggplot2::scale_color_gradient(
        low = "lightgray",
        high = "purple"),
    colMin = -2.5,
    colMax = 2.5,
    dotMin = 0L,
    dotScale = 6L) {
    .checkFormat(format)
    if (format == "ensgene") {
        genes <- .convertGenesToSymbols(object, genes = genes)
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
            !!genes) %>%
        group_by(!!!syms(c("ident", "gene"))) %>%
        summarize(
            avgExp = mean(expm1(.data[["expression"]])),
            pctExp = .percentAbove(.data[["expression"]], threshold = 0L)
        ) %>%
        ungroup() %>%
        group_by(!!sym("gene")) %>%
        mutate(
            avgExpScale = scale(.data[["avgExp"]]),
            avgExpScale = .minMax(
                .data[["avgExpScale"]],
                max = colMax,
                min = colMin)
        )
    data[["pctExp"]][data[["pctExp"]] < dotMin] <- NA

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "gene",
            y = "ident")
    ) +
        geom_point(
            mapping = aes_string(
                color = "avgExpScale",
                size = "pctExp")
        ) +
        scale_radius(range = c(0L, dotScale)) +
        labs(x = "gene",
             y = "ident")
    if (is(color, "ScaleContinuous")) {
        p <- p + color
    }

    p
}




# Methods ======================================================================
#' @rdname plotDot
#' @export
setMethod(
    "plotDot",
    signature("seurat"),
    .plotDot.seurat)
