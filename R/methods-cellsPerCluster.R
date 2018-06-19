#' Cell Counts per Cluster
#'
#' @name cellsPerCluster
#'
#' @inheritParams general
#' @param return Return as `list` or `markdown`.
#'
#' @examples
#' # seurat ====
#' # list
#' x <- cellsPerCluster(seurat_small, return = "list")
#' names(x)
#'
#' # markdown
#' cellsPerCluster(seurat_small, return = "markdown")
NULL



# Constructors =================================================================
.cellsPerClusterMarkdown <- function(list, level = 2L) {
    invisible(mapply(
        x = list,
        title = names(list),
        FUN = function(x, title) {
            markdownHeader(title, level = level, asis = TRUE)
            show(kable(
                x = x,
                caption = title,
                digits = 2L
            ))
        }
    ))
}



# Methods ======================================================================
#' @rdname cellsPerCluster
#' @export
setMethod(
    "cellsPerCluster",
    signature("seurat"),
    function(
        object,
        interestingGroups,
        return = c("list", "markdown")
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        return <- match.arg(return)

        ident <- slot(object, "ident")
        assert_is_factor(ident)

        list <- mclapply(
            X = levels(ident),
            FUN = function(ident) {
                subset <- SubsetData(object, ident.use = ident)
                metrics <- metrics(
                    object = subset,
                    interestingGroups = interestingGroups
                )
                assert_is_subset("interestingGroups", colnames(metrics))
                metrics %>%
                    arrange(!!!syms(interestingGroups)) %>%
                    group_by(!!!syms(interestingGroups)) %>%
                    summarize(n = n()) %>%
                    mutate(ratio = !!sym("n") / sum(!!sym("n"))) %>%
                    group_by(!!!syms(interestingGroups))
            }
        )
        names(list) <- levels(ident)

        if (return == "list") {
            list
        } else if (return == "markdown") {
            .cellsPerClusterMarkdown(list)
        }
    }
)
