#' Create a new Seurat object
#'
#' @author Michael Steinbaugh
#'
#' @param object [SummarizedExperiment] containing filtered counts.
#'
#' @return [seurat].
#' @export
new_seurat <- function(object, ...) {
    # Filtering criteria
    min.genes <- metadata(data) %>%
        .[["filtering_criteria"]] %>%
        .[["genes"]]
    max.mito.ratio <- metadata(data) %>%
        .[["filtering_criteria"]] %>%
        .[["mito_ratio"]]

    seurat <- new("seurat", raw.data = counts(data, gene2symbol = TRUE)) %>%
        Setup(meta.data = metrics(data) %>% dotted(rownames = FALSE),
              min.genes = min.genes,
              ...) %>%
        SubsetData(subset.name = "mito.ratio", accept.high = max.mito.ratio)

    message(paste("Seurat object:", object_size(seurat)))
    seurat
}
