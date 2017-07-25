#' S4 Generics
#'
#' @rdname all_generics
#' @name all_generics
#' @keywords internal
#'
#' @param object Object.
#' @param x Primary object.
#' @param y Secondary object.
#' @param ... Additional arguments.
NULL



#' @rdname aggregate_replicates
#' @inherit all_generics
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



#' @rdname bcbio
#' @inherit all_generics
#' @export
setGeneric("bcbio", function(object, ...) {
    standardGeneric("bcbio")
})



#' @rdname bcbio
#' @usage NULL
#' @export
setGeneric("bcbio<-", function(object, ..., value) {
    standardGeneric("bcbio<-")
})



#' @rdname metadata_table
#' @inherit all_generics
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})



#' @rdname interesting_groups
#' @inherit all_generics
#' @export
setGeneric("interesting_groups", function(object) {
    standardGeneric("interesting_groups")
})



#' @rdname known_markers_detected
#' @inherit all_generics
#' @export
setGeneric("known_markers_detected", function(x, y, ...) {
    standardGeneric("known_markers_detected")
})



#' @rdname metrics
#' @inherit all_generics
#' @export
setGeneric("metrics", function(object, ...) {
    standardGeneric("metrics")
})



#' @rdname plot_cell_counts
#' @inherit all_generics
#' @export
setGeneric("plot_cell_counts", function(object, ...) {
    standardGeneric("plot_cell_counts")
})



#' @rdname plot_cellular_barcodes
#' @inherit all_generics
#' @export
setGeneric("plot_cellular_barcodes", function(object, ...) {
    standardGeneric("plot_cellular_barcodes")
})



#' @rdname plot_clusters
#' @inherit all_generics
#' @export
setGeneric("plot_clusters", function(object, ...) {
    standardGeneric("plot_clusters")
})



#' @rdname plot_genes_per_cell
#' @inherit all_generics
#' @export
setGeneric("plot_genes_per_cell", function(object, ...) {
    standardGeneric("plot_genes_per_cell")
})



#' @rdname plot_known_markers
#' @inherit all_generics
#' @export
setGeneric("plot_known_markers", function(x, y, ...) {
    standardGeneric("plot_known_markers")
})



#' @rdname plot_mito_ratio
#' @inherit all_generics
#' @export
setGeneric("plot_mito_ratio", function(object, ...) {
    standardGeneric("plot_mito_ratio")
})



#' @rdname plot_novelty
#' @inherit all_generics
#' @export
setGeneric("plot_novelty", function(object, ...) {
    standardGeneric("plot_novelty")
})



#' @rdname plot_top_markers
#' @inherit all_generics
#' @export
setGeneric("plot_top_markers", function(x, y, ...) {
    standardGeneric("plot_top_markers")
})



#' @rdname plot_umis_per_cell
#' @inherit all_generics
#' @export
setGeneric("plot_umis_per_cell", function(object, ...) {
    standardGeneric("plot_umis_per_cell")
})



#' @rdname plot_umis_vs_genes
#' @inherit all_generics
#' @export
setGeneric("plot_umis_vs_genes", function(object, ...) {
    standardGeneric("plot_umis_vs_genes")
})



#' @rdname sample_metadata
#' @inherit all_generics
#' @export
setGeneric("sample_metadata", function(object, ...) {
    standardGeneric("sample_metadata")
})



#' @rdname select_samples
#' @inherit all_generics
#' @export
setGeneric("select_samples", function(object, ...) {
    standardGeneric("select_samples")
})



#' @rdname top_barcodes
#' @inherit all_generics
#' @export
setGeneric("top_barcodes", function(object, ...) {
    standardGeneric("top_barcodes")
})



#' @rdname top_markers
#' @inherit all_generics
#' @export
setGeneric("top_markers", function(object, ...) {
    standardGeneric("top_markers")
})
