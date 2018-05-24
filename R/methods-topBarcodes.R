#' Top Barcodes
#'
#' @name topBarcodes
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param n Number of barcodes to return per sample.
#'
#' @return `grouped_df` containing [metrics()] return. Grouped by `sampleID`
#'   and arranged by `nUMI` column in descending order.
#'
#' @examples
#' # SingleCellExperiment ====
#' x <- topBarcodes(cellranger_small)
#' glimpse(x)
NULL



# Methods ======================================================================
#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("SingleCellExperiment"),
    function(object, n = 10L) {
        validObject(object)
        metrics <- metrics(object)
        assert_is_subset("nUMI", colnames(metrics))
        metrics %>%
            rownames_to_column("cellID") %>%
            arrange(desc(!!sym("nUMI"))) %>%
            group_by(!!sym("sampleID")) %>%
            slice(seq_len(n))
    }
)




#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("seurat"),
    getMethod("topBarcodes", "SingleCellExperiment")
)
