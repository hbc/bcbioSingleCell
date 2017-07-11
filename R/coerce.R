#' Object coercion definitions
#'
#' @rdname coerce
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @param from Object to coerce.



setAs("SCSubset", "seurat", function(from) {
    # Filtering criteria
    min_genes <- metadata(from) %>%
        .[["filtering_criteria"]] %>%
        .[["min_genes"]]
    max_mito_ratio <- metadata(from) %>%
        .[["filtering_criteria"]] %>%
        .[["max_mito_ratio"]]
    new("seurat", raw.data = counts(from, gene2symbol = TRUE)) %>%
        Setup(meta.data = metrics(from) %>% dotted(rownames = FALSE),
              min.genes = min_genes,
              # Use the object name as the Seurat project name
              project = deparse(substitute(from))) %>%
        SubsetData(subset.name = "mito.ratio", accept.high = max_mito_ratio)
})
