#' Fetch t-SNE Locations and Cellular Metadata
#'
#' @rdname fetchTSNEData
#' @name fetchTSNEData
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @return [data.frame] of t-SNE points and metadata for each cell.
#'
#' @examples
#' \dontrun{
#' data(seurat)
#' tsne <- fetchTSNEData(seurat)
#' }
NULL



# Methods ====
#' @rdname fetchTSNEData
#' @export
setMethod("fetchTSNEData", "seurat", function(object) {
    dimCode <- c("tSNE_1", "tSNE_2")
    FetchData(object, dimCode) %>%
        rownames_to_column("cell") %>%
        left_join(object@meta.data %>%
                      rownames_to_column("cell"),
                  by = "cell") %>%
        left_join(data.frame(object@ident) %>%
                      rownames_to_column("cell"),
                  by = "cell") %>%
        group_by(!!sym("object.ident")) %>%
        mutate(centerx = median(.data[["tSNE_1"]]),
               centery = median(.data[["tSNE_2"]]))
})
