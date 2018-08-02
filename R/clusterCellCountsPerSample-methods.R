#' Cluster Cell Counts per Sample
#'
#' @name clusterCellCountsPerSample
#'
#' @inheritParams general
#'
#' @return `tibble` grouped by "`sampleName" column, arranged by abundance.
#'
#' @examples
#' # SingleCellExperiment ====
#' clusterCellCountsPerSample(indrops_small)
NULL



# Methods ======================================================================
#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    "clusterCellCountsPerSample",
    signature("SingleCellExperiment"),
    function(object) {
        .assertHasIdent(object)
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
