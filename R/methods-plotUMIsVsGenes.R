#' Plot UMI and Gene Correlation
#'
#' @rdname plotUMIsVsGenes
#' @name plotUMIsVsGenes
#'
#' @return [ggplot].
NULL



# Constructors ====
.plotUMIsVsGenes <- function(object) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~nUMI,
             y = ~nGene,
             color = ~sampleName)) +
        labs(x = "umis per cell",
             y = "genes per cell") +
        geom_point(size=0.8, alpha=0.3) +
        geom_smooth(se=FALSE) +
        scale_x_log10() +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



# Methods ====
#' @rdname plotUMIsVsGenes
#' @export
setMethod("plotUMIsVsGenes", "bcbioSCDataSet", .plotUMIsVsGenes)



#' @rdname plotUMIsVsGenes
#' @export
setMethod("plotUMIsVsGenes", "bcbioSCSubset", .plotUMIsVsGenes)
