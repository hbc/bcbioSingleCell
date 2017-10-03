#' Fetch Data
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param object [seurat] object.
#' @param dimCode Character vector of X and Y coordinate data to be used for
#'   plotting. This can be `c("tSNE_1", "tSNE_2")` for tSNE data, or `c("PC1",
#'   "PC2")` for PCA data.
#'
#' @return [data.frame].
.fetchDimDataSeurat <- function(object, dimCode) {
    meta <- object@meta.data %>%
        camel(strict = FALSE) %>%
        rownames_to_column("cell")
    ident <- object@ident %>%
        as.data.frame() %>%
        set_colnames("ident") %>%
        rownames_to_column("cell")
    FetchData(object, vars.all = dimCode) %>%
        as.data.frame() %>%
        camel(strict = FALSE) %>%
        rownames_to_column("cell") %>%
        left_join(meta, by = "cell") %>%
        left_join(ident, by = "cell") %>%
        group_by(!!sym("ident")) %>%
        mutate(centerX = median(.data[[camel(dimCode[[1]], strict = FALSE)]]),
               centerY = median(.data[[camel(dimCode[[2]], strict = FALSE)]]))
}
