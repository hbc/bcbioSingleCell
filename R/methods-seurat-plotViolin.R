# TODO Add facet return method



#' Plot Violin
#'
#' @name plotViolin
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotViolin
#'
#' @inheritParams general
#' @inheritParams plotDot
#' @inheritParams ggplot2::geom_violin
#'
#' @examples
#' # seurat ====
#' genes <- head(rownames(Seurat::pbmc_small), 2L)
#' print(genes)
#'
#' # grid
#' plotViolin(Seurat::pbmc_small, genes = genes, return = "grid")
#'
#' # list
#' list <- plotViolin(Seurat::pbmc_small, genes = genes, return = "list")
#' names(list)
#'
#' # markdown
#' plotViolin(Seurat::pbmc_small, genes = genes, return = "markdown")
NULL



# Methods ======================================================================
#' @rdname plotViolin
#' @export
setMethod(
    "plotViolin",
    signature("seurat"),
    function(
        object,
        genes,
        scale = c("count", "width", "area"),
        fill = scale_fill_hue(),
        dark = FALSE,
        headerLevel = 2L,
        return = c("grid", "list", "markdown")
    ) {
        scale <- match.arg(scale)
        assert_is_any_of(fill, c("ScaleDiscrete", "character", "NULL"))
        if (is.character(fill)) {
            assert_is_a_string(fill)
        }
        assert_is_a_bool(dark)
        assertIsAHeaderLevel(headerLevel)
        return <- match.arg(return)

        # Fetch data ===========================================================
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

        # Loop across genes ====================================================
        plotlist <- lapply(seq_along(genes), function(a) {
            gene <- genes[[a]]
            data <- data[data[["gene"]] == gene, , drop = FALSE]

            # Dynamically provide mapped or single color support
            if (is_a_string(fill)) {
                fillArg <- fill
            } else {
                fillArg <- NULL
            }

            if (isTRUE(dark)) {
                color <- NA
            } else {
                color <- "black"
            }

            violin <- geom_violin(
                mapping = aes_string(fill = "ident"),
                # never include a color border
                color = color,
                scale = scale,
                adjust = 1L,
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
                labs(title = gene) +
                guides(fill = FALSE)

            if (is(fill, "ScaleDiscrete")) {
                p <- p + fill
            }

            if (isTRUE(dark)) {
                p <- p + theme_midnight()
            } else {
                p <- p + theme_paperwhite()
            }

            p
        })
        names(plotlist) <- genes

        # Return ===============================================================
        dynamicPlotlist(
            plotlist,
            return = return,
            headerLevel = headerLevel
        )
    }
)
