##' Fetch t-SNE locations and cellular metadata
##'
##' @param seurat a seurat object
##' @return dataframe of t-SNE points and metadata for each cell
##' @author Rory Kirchner
FetchTsneData = function(seurat) {
  dim.code = c("tSNE_1", "tSNE_2")
  FetchData(seurat, dim.code) %>%
    tibble::rownames_to_column("cell") %>%
    left_join(seurat@meta.data %>%
              tibble::rownames_to_column("cell")) %>%
    left_join(data.frame(seurat@ident) %>%
              tibble::rownames_to_column("cell")) %>%
    group_by(seurat.ident) %>%
    mutate(centerx = median(tSNE_1), centery = median(tSNE_2))
}

##' Fetch t-SNE locations, cellular metadata and expression of genes
##'
##' This gets t-SNE locations, cellular metadata, expression of genes and
##' the geometric mean of the gene expression from a seurat object.
##'
##' @param seurat a seurat object
##' @param genes genes of which to get expression data
##' @return tidy dataframe of t-SNE points, cellular metadata and gene expression
##' @author Rory Kirchner
FetchTsneExpressionData = function(seurat, genes) {
  dat = FetchData(seurat, genes) %>%
    as.data.frame()
  dat$geomean = colMeans(t(dat))
  dat = dat %>%
    tibble::rownames_to_column("cell")
  FetchTsneData(seurat) %>%
    left_join(dat) %>%
    tidyr::gather(gene, expresssion, !!genes)
}

##' Plot fetched t-SNE expression data
##'
##' @param df a dataframe from FetchTsneExpressionData
##' @return ggplot2 t-SNE plot
##' @author Rory Kirchner
PlotTsneExpressionData = function(df) {
  ggplot(df, aes(tSNE_1, tSNE_2, color=geomean)) +
    geom_text(aes(tSNE_1, tSNE_2, label=seurat.ident, color=geomean), size=2) +
    scale_color_viridis() +
    guides(color=FALSE) +
    geom_text(aes(x=centerx, y=centery, label=seurat.ident), color="red") +
    xlab("") +
    ylab("") +
    DarkTheme() +
    theme(text=element_text(size=8))
}
