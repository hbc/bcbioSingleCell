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
#'
#' @examples
#' # seurat ====
#' genes <- slot(pbmc_small, "raw.data") %>% rownames() %>% head(2L)
#'
#' # grid
#' plotViolin(pbmc_small, genes = genes, return = "grid")
#'
#' # list
#' list <- plotViolin(pbmc_small, genes = genes, return = "list")
#' names(list)
#'
#' # markdown
#' plotViolin(pbmc_small, genes = genes, return = "markdown")
NULL



# Constructors =================================================================
.plotViolin.seurat <- function(  # nolint
    object,
    genes,
    fill = scale_fill_viridis(discrete = TRUE),
    headerLevel = 2L,
    return = c("grid", "list", "markdown")
) {
    assertIsFillScaleDiscreteOrNULL(fill)
    assertIsAHeaderLevel(headerLevel)
    return <- match.arg(return)

    # Fetch data ===============================================================
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
        p <- ggplot(
            data,
            mapping = aes_string(
                x = "ident",
                y = "expression",
                fill = "ident"
            )
        ) +
            geom_violin(
                color = "black",
                scale = "width",
                adjust = 1L,
                trim = TRUE
            ) +
            labs(title = gene) +
            guides(fill = FALSE)
        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
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



# Methods ======================================================================
#' @rdname plotViolin
#' @export
setMethod(
    "plotViolin",
    signature("seurat"),
    .plotViolin.seurat
)
