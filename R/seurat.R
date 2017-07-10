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
        .[["min_genes"]]
    max_mito_ratio <- metadata(object) %>%
        .[["filtering_criteria"]] %>%
        .[["max_mito_ratio"]]
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
#' @param seurat [seurat].
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



#' Plot known markers
#'
#' @param seurat [seurat].
#' @param cell_type Cell type.
#'
#' @export
plot_known_markers <- function(seurat, cell_type) {
    lapply(seq_along(cell_type), function(a) {
        writeLines(
            c("", "",
              paste("###", cell_type[a]),
              "", ""))
        genes <- known_markers_detected %>%
            .[.$cell_type == cell_type[a], "gene"] %>% unique %>% sort
        plot_clusters(seurat, genes)
    }) %>% invisible
}



#' Plot top markers
#'
#' @param seurat [seurat].
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
