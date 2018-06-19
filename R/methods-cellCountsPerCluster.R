#' Cell Counts per Cluster
#'
#' @name cellCountsPerCluster
#'
#' @inheritParams general
#'
#' @return `tibble` grouped by "`ident`" column, arranged by abundance.
#'
#' @examples
#' # seurat ====
#' # list
#' x <- cellCountsPerCluster(seurat_small, return = "list")
#' names(x)
#'
#' # markdown
#' cellCountsPerCluster(seurat_small, return = "markdown")
NULL



# Methods ======================================================================
#' @rdname cellCountsPerCluster
#' @export
setMethod(
    "cellCountsPerCluster",
    signature("seurat"),
    function(object, interestingGroups) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        cols <- unique(c("ident", interestingGroups))
        assert_is_subset(cols, colnames(metrics))
        metrics %>%
            arrange(!!!syms(cols)) %>%
            group_by(!!!syms(cols)) %>%
            summarize(n = n()) %>%
            arrange(desc(!!sym("n")), .by_group = TRUE) %>%
            mutate(ratio = !!sym("n") / sum(!!sym("n")))
    }
)
