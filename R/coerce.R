#' Object coercion definitions
#'
#' @rdname coerce
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @param from Object to coerce.



setAs("SCSubset", "seurat", function(from) {
    name <- deparse(substitute(from))
    raw_data <- counts(from, gene2symbol = TRUE)
    meta_data <- metrics(from) %>% dotted(rownames = FALSE)

    # Filtering criteria
    min_genes <- metadata(from) %>%
        .[["filtering_criteria"]] %>%
        .[["min_genes"]]
    max_mito_ratio <- metadata(from) %>%
        .[["filtering_criteria"]] %>%
        .[["max_mito_ratio"]]

    new("seurat", raw.data = raw_data) %>%
        Setup(meta.data = meta_data,
              min.genes = min_genes,
              project = name) %>%
        SubsetData(subset.name = "mito.ratio", accept.high = max_mito_ratio)
})
