#' Create a new Seurat object
#'
#' @author Michael Steinbaugh
#'
#' @param object [SummarizedExperiment] containing filtered counts.
#' @param ... Additional parameters, passed to [Seurat::Setup()].
#'
#' @return [seurat].
#' @export
new_seurat <- function(object, ...) {
    # Filtering criteria
    min_genes <- metadata(object) %>%
        .[["filtering_criteria"]] %>%
        .[["genes"]]
    max_mito_ratio <- metadata(object) %>%
        .[["filtering_criteria"]] %>%
        .[["mito_ratio"]]
    seurat <- new("seurat", raw.data = counts(object, gene2symbol = TRUE)) %>%
        Setup(meta.data = metrics(object) %>% dotted(rownames = FALSE),
              min.genes = min_genes,
              ...) %>%
        SubsetData(subset.name = "mito.ratio", accept.high = max_mito_ratio)
    message(paste("Seurat object:", object_size(seurat)))
    seurat
}



#' Plot clusters
#'
#' @param seurat Seurat object.
#' @param genes Genes.
#'
#' @return [ggplot].
#' @export
plot_clusters <- function(seurat, genes) {
    VlnPlot(seurat, genes,
            use.scaled = TRUE)
    FeaturePlot(seurat, genes,
                reduction.use = "tsne",
                cols.use = c("grey", "blue"))
    FeaturePlot(seurat, genes,
                reduction.use = "pca",
                cols.use = c("grey", "red"))
}



#' Top marker genes
#'
#' @param markers Markers.
#' @param n Number of genes per cluster.
#'
#' @return [data.frame]
#' @export
top_markers <- function(markers, n = 4L) {
    markers %>%
        group_by(.data[["cluster"]]) %>%
        # FIXME Need SE version heres
        top_n(n = n, wt = .data[["avg_diff"]])
}



#' Plot top markers
#'
#' @param seurat Seurat object.
#' @param top Top genes data frame.
#'
#' @return [ggplot].
#' @export
plot_top_markers <- function(seurat, top) {
    lapply(levels(top[["cluster"]]), function(a) {
        writeLines(
            c("", "",
              paste("###", "Cluster", a),
              "", ""))
        genes <- top %>%
            .[.[["cluster"]] == a, "gene"]
        plot_clusters(seurat, genes)
    }) %>% invisible
}
