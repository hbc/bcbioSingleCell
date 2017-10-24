#' Plot Features
#'
#' @rdname plotFeatures
#' @name plotFeatures
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param features Character vector of features (e.g. gene expression, PC
#'   scores, number of genes detected).
#'
#' @return No return, only graphical output.
NULL



# Methods ====
#' @rdname plotFeatures
#' @importFrom Seurat FeaturePlot
#' @export
setMethod("plotFeatures", "seurat", function(object, features) {
    FeaturePlot(
        object,
        features.plot = features,
        cols.use = c("white", "black"),
        no.legend = TRUE,
        do.return = FALSE
    )
})
