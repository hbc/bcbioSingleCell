#' Plot Violin
#'
#' @name plotViolin
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotViolin
#'
#' @inheritParams general
#' @inheritParams plotDot
#'
#' @param fill Fill color palette. Defaults to viridis.
#' @param return Return type. "grid", "list", and "markdown" are supported.
#' @param headerLevel Markdown header level. Only applicable when
#'   `return = "markdown"`.
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
#' @importFrom basejump dynamicPlotlist
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string geom_violin ggplot
#' @importFrom rlang !! sym
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom viridis scale_fill_viridis
.plotViolin.seurat <- function(  # nolint
    object,
    genes,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    return = "grid",
    headerLevel = 2L
) {
    validReturn <- c("grid", "list", "markdown")
    if (!return %in% validReturn) {
        abort(paste("`return` must contain:", toString(validReturn)))
    }

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
