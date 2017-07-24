#' Object coercion definitions
#'
#' @rdname coerce
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @param from Object to coerce.



# https://goo.gl/GMQAdC
setAs("SCSubset", "seurat", function(from) {
    name <- deparse(substitute(from))
    counts <- counts(from, gene2symbol = TRUE)
    metrics <- metrics(from) %>% dotted

    # Filtering criteria
    min_genes <- metadata(from) %>%
        .[["filtering_criteria"]] %>%
        .[["min_genes"]]
    max_mito_ratio <- metadata(from) %>%
        .[["filtering_criteria"]] %>%
        .[["max_mito_ratio"]]

    new("seurat", raw.data = counts) %>%
        Setup(project = name,
              min.genes = min_genes) %>%
        AddMetaData(metrics) %>%
        SubsetData(subset.name = "mito.ratio", accept.high = max_mito_ratio)
})
