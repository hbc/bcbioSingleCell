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
#' object <- indrops_small
#'
#' # tibble
#' x <- topBarcodes(object, return = "tibble")
#' glimpse(x)
#'
#' # list
#' x <- topBarcodes(object, return = "list")
#' lapply(x, class)
NULL



.topBarcodes.SCE <-  # nolint
    function(
        object,
        n = 100L,
        return = c("tibble", "list")
    ) {
        validObject(object)
        assertIsAnImplicitInteger(n)
        return <- match.arg(return)

        metrics <- metrics(object)
        cols <- c("cellID", "sampleID", "sampleName", "nUMI")
        assert_is_subset(cols, colnames(metrics))
        data <- metrics %>%
            select(!!!syms(cols)) %>%
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



#' @rdname topBarcodes
#' @export
setMethod(
    f = "topBarcodes",
    signature = signature("SingleCellExperiment"),
    definition = .topBarcodes.SCE
)
