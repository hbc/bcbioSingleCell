#' S4 Generics
#'
#' @rdname all_generics
#' @name all_generics
#' @keywords internal
#'
#' @param object Object.
#' @param x Primary object.
#' @param y Secondary object.
#' @param value Object to assign.
#' @param ... Additional arguments.
#'
#' @return No value.
NULL



#' @rdname aggregate_replicates
#' @inheritParams all_generics
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



#' @rdname bcbio
#' @inheritParams all_generics
#' @export
setGeneric("bcbio", function(object, ...) {
    standardGeneric("bcbio")
})



#' @rdname bcbio
#' @inheritParams all_generics
#' @export
setGeneric("bcbio<-", function(object, ..., value) {
    standardGeneric("bcbio<-")
})



#' @rdname filter_cells
#' @inheritParams all_generics
#' @export
setGeneric("filter_cells", function(object, ...) {
    standardGeneric("filter_cells")
})



#' @rdname metadata_table
#' @inheritParams all_generics
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})



#' @rdname interesting_groups
#' @inheritParams all_generics
#' @export
setGeneric("interesting_groups", function(object) {
    standardGeneric("interesting_groups")
})



#' @rdname known_markers_detected
#' @inheritParams all_generics
#' @export
setGeneric("known_markers_detected", function(x, y, ...) {
    standardGeneric("known_markers_detected")
})



#' @rdname metrics
#' @inheritParams all_generics
#' @export
setGeneric("metrics", function(object, ...) {
    standardGeneric("metrics")
})



#' @rdname plot_cell_counts
#' @inheritParams all_generics
#' @export
setGeneric("plot_cell_counts", function(object, ...) {
    standardGeneric("plot_cell_counts")
})



#' @rdname plot_cellular_barcodes
#' @inheritParams all_generics
#' @export
setGeneric("plot_cellular_barcodes", function(object, ...) {
    standardGeneric("plot_cellular_barcodes")
})



#' @rdname plot_clusters
#' @inheritParams all_generics
#' @export
setGeneric("plot_clusters", function(object, ...) {
    standardGeneric("plot_clusters")
})



#' @rdname plot_genes_per_cell
#' @inheritParams all_generics
#' @export
setGeneric("plot_genes_per_cell", function(object, ...) {
    standardGeneric("plot_genes_per_cell")
})



#' @rdname plot_known_markers
#' @inheritParams all_generics
#' @export
setGeneric("plot_known_markers", function(x, y, ...) {
    standardGeneric("plot_known_markers")
})



#' @rdname plot_mito_ratio
#' @inheritParams all_generics
#' @export
setGeneric("plot_mito_ratio", function(object, ...) {
    standardGeneric("plot_mito_ratio")
})



#' @rdname plot_novelty
#' @inheritParams all_generics
#' @export
setGeneric("plot_novelty", function(object, ...) {
    standardGeneric("plot_novelty")
})



#' @rdname plot_top_markers
#' @inheritParams all_generics
#' @export
setGeneric("plot_top_markers", function(x, y, ...) {
    standardGeneric("plot_top_markers")
})



#' @rdname plot_umis_per_cell
#' @inheritParams all_generics
#' @export
setGeneric("plot_umis_per_cell", function(object, ...) {
    standardGeneric("plot_umis_per_cell")
})



#' @rdname plot_umis_vs_genes
#' @inheritParams all_generics
#' @export
setGeneric("plot_umis_vs_genes", function(object, ...) {
    standardGeneric("plot_umis_vs_genes")
})



#' @rdname sample_metadata
#' @inheritParams all_generics
#' @export
setGeneric("sample_metadata", function(object, ...) {
    standardGeneric("sample_metadata")
})



#' @rdname select_samples
#' @inheritParams all_generics
#' @export
setGeneric("select_samples", function(object, ...) {
    standardGeneric("select_samples")
})



#' @rdname top_barcodes
#' @inheritParams all_generics
#' @export
setGeneric("top_barcodes", function(object, ...) {
    standardGeneric("top_barcodes")
})



#' @rdname top_markers
#' @inheritParams all_generics
#' @export
setGeneric("top_markers", function(object, ...) {
    standardGeneric("top_markers")
})
