##' Plot fetched t-SNE expression data
##'
##' @param df a dataframe from FetchTsneExpressionData
##' @param colorpoints color points by geometric mean or expression of
##'  individual gene
##' @return ggplot2 t-SNE plot
##' @examples
##' \dontrun{
##'   dat = FetchTsneExpressionData(seurat, topGenes)
##'   # to make t-SNE plot colored by geometric mean of topGenes
##'   PlotTsneExpressionData(dat, colorpoints="geometricmean")
##'   # to make faceted t-SNE plot of each gene (looks good at up to 6 genes)
##'   PlotTsneExpression(dat, colorpoints="expression") +
##'     facet_wrap(~gene)
##' } @author Rory Kirchner
PlotTsneExpressionData = function(df, colorpoints="geomean") {
  ggplot(df, aes_string("tSNE_1", "tSNE_2", color=colorpoints)) +
    geom_text(aes_string("tSNE_1", "tSNE_2", label="seurat.ident",
                         color=colorpoints),
              size=2) +
    scale_color_viridis() +
    guides(color=FALSE) +
    geom_text(aes(x=centerx, y=centery, label=seurat.ident), color="pink", size=3,
              fontface="bold") +
    xlab("") +
    ylab("") +
    DarkTheme() +
    theme(strip.background=element_rect(fill="black"),
          strip.text = element_text(face="bold"),
          text=element_text(size=8))
}
