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
#' - "`tibble`": `grouped_df`. Grouped by `sampleID` and arranged by `nUMI`
#'   column in descending order. Cellular barcodes are in the `cellID` column.
#' - "`list`": `list`. Top barcodes as `character`, split by `sampleID`.
#'
#' @examples
#' # tibble
#' x <- topBarcodes(indrops_small, return = "tibble")
#' glimpse(x)
#'
#' # list
#' x <- topBarcodes(indrops_small, return = "list")
#' lapply(x, class)
NULL



#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("SingleCellExperiment"),
    function(
        object,
        n = 100L,
        return = c("tibble", "list")
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

        if (return == "tibble") {
            data
        } else if (return == "list") {
            data %>%
                split(.[["sampleID"]]) %>%
                map("cellID")
        }
    }
)
