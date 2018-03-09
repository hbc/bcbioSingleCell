#' Add Median Labels to a Plot
#'
#' For example, `medianCol` can be `nGene`. Median values are always calculated
#' by sample (`sampleName`).
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom ggplot2 aes_string geom_label
#' @importFrom S4Vectors aggregate
#'
#' @return ggplot2 `geom_label`.
.medianLabels <- function(metrics, medianCol, digits = 0L) {
    formula <- formula(paste(medianCol, "sampleName", sep = " ~ "))
    data <- aggregate(
        formula = formula,
        data = metrics,
        FUN = median
    )
    data[["roundedMedian"]] <- round(data[[medianCol]], digits = digits)

    # Add `sampleNameAggregate` column for facet wrapping, if necessary
    if ("sampleNameAggregate" %in% colnames(metrics)) {
        sampleFacet <- metrics[, c("sampleName", "sampleNameAggregate")] %>%
        distinct()
        data <- left_join(data, sampleFacet, by = "sampleName")
    }

    geom_label(
        data = data,
        mapping = aes_string(label = "roundedMedian"),
        alpha = qcLabelAlpha,
        color = qcLabelColor,
        fill = qcLabelFill,
        fontface = qcLabelFontface,
        label.padding = qcLabelPadding,
        label.size = qcLabelSize,
        show.legend = FALSE
    )
}
