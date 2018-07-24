#' Top Barcodes
#'
#' Obtain the top cellular barcodes based on UMI counts.
#'
#' @name topBarcodes
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param n `scalar integer`. Number of barcodes to return per sample.
#'
#' @return
#' - "`list`": `list` containing top barcodes as a `character` vector, split by
#'   `sampleID`.
#' - "`data.frame`": `tibble` grouped by `sampleID` and arranged by `nUMI`
#'   column in descending order. Cellular barcodes are in the `cellID` column.
#'
#' @examples
#' # SingleCellExperiment ====
#' # list
#' x <- topBarcodes(cellranger_small, return = "list")
#' lapply(x, class)
#'
#' # data.frame
#' x <- topBarcodes(cellranger_small, return = "data.frame")
#' glimpse(x)
NULL



# Methods ======================================================================
#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("SingleCellExperiment"),
    function(
        object,
        n = 1000L,
        return = c("data.frame", "list")
    ) {
        validObject(object)
        assertIsAnImplicitInteger(n)
        return <- match.arg(return)

        metrics <- metrics(object)
        cols <- c("sampleID", "sampleName", "nUMI")
        assert_is_subset(cols, colnames(metrics))
        data <- metrics %>%
            as.data.frame() %>%
            rownames_to_column("cellID") %>%
            select(!!!syms(c(cols, "cellID"))) %>%
            group_by(!!sym("sampleID")) %>%
            arrange(desc(!!sym("nUMI")), .by_group = TRUE) %>%
            slice(seq_len(n))

        if (return == "data.frame") {
            data
        } else if (return == "list") {
            data %>%
                split(.[["sampleID"]]) %>%
                map("cellID")
        }
    }
)



#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("seurat"),
    getMethod("topBarcodes", "SingleCellExperiment")
)
