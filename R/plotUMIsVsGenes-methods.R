#' Plot UMI and Gene Correlation
#'
#' @rdname plotUMIsVsGenes
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit plotGenesPerCell
#'
#' @return [ggplot].



#' @rdname plotUMIsVsGenes
#' @usage NULL
.plotUMIsVsGenes <- function(object) {
    metrics <- metrics(object)
    p <- ggplot(metrics,
        aes_(x = ~umiCounts,
             y = ~genesDetected,
             color = ~sampleName)) +
        labs(x = "umis per cell",
             y = "genes per cell") +
        geom_smooth(method = "lm", se = FALSE) +
        scale_x_log10() +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



#' @rdname plotUMIsVsGenes
#' @export
setMethod("plotUMIsVsGenes", "bcbioSCDataSet", .plotUMIsVsGenes)



#' @rdname plotUMIsVsGenes
#' @export
setMethod("plotUMIsVsGenes", "bcbioSCSubset", .plotUMIsVsGenes)
