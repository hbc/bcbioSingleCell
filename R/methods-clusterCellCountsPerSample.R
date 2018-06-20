#' Cluster Cell Counts per Sample
#'
#' @name clusterCellCountsPerSample
#'
#' @inheritParams general
#'
#' @return `tibble` grouped by "`sampleName" column, arranged by abundance.
#'
#' @examples
#' # seurat ====
#' clusterCellCountsPerSample(seurat_small)
NULL



# Methods ======================================================================
#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    "clusterCellCountsPerSample",
    signature("SingleCellExperiment"),
    function(object) {
        metrics <- metrics(object)
        cols <- c("sampleName", "ident")
        assert_is_subset(cols, colnames(metrics))
        metrics %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!!syms(cols)) %>%
            summarize(n = n()) %>%
            ungroup() %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!sym("sampleName")) %>%
            mutate(ratio = !!sym("n") / sum(!!sym("n")))
    }
)



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    "clusterCellCountsPerSample",
    signature("seurat"),
    getMethod("clusterCellCountsPerSample", "SingleCellExperiment")
)
