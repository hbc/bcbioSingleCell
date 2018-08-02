#' Plot Violin
#'
#' @name plotViolin
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotViolin
#'
#' @inheritParams general
#' @inheritParams ggplot2::geom_violin
#'
#' @seealso [Seurat::VlnPlot()].
#'
#' @return `ggplot`.
#'
#' @examples
#' # SingleCellExperiment ====
#' object <- cellranger_small
#' genes <- head(rownames(object), n = 4L)
#' glimpse(genes)
#' plotViolin(object, genes = genes)
NULL



# Methods ======================================================================
#' @rdname plotViolin
#' @export
setMethod(
    "plotViolin",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        scale = c("count", "width", "area"),
        fill = getOption("bcbio.discrete.fill", NULL),
        legend = getOption("bcbio.legend", TRUE)
    ) {
        scale <- match.arg(scale)
        assert_is_any_of(fill, c("ScaleDiscrete", "character", "NULL"))
        if (is.character(fill)) {
            assert_is_a_string(fill)
        }

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
            group_by(!!sym("gene"))

        violin <- geom_violin(
            mapping = aes(fill = !!sym("ident")),
            # never include a color border
            color = "black",
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
            mapping = aes(
                x = !!sym("ident"),
                y = !!sym("expression")
            )
        ) +
            violin +
            # Note that `scales = free_y` will hide the x-axis for some plots
            facet_wrap(facets = sym("gene"), scales = "free_y")

        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }

        p
    }
)
