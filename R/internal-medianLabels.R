#' Add Median Labels to a Plot
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom ggplot2 aes_string geom_label
#' @importFrom S4Vectors aggregate
.medianLabels <- function(metrics, medianCol, digits = 0) {
    # For example, `medianCol` can be `nGene`.
    # Median values are always calculated by sample (`sampleName`).
    formula <- formula(paste(medianCol, "sampleName", sep = " ~ "))
    data <- aggregate(
        formula = formula,
        data = metrics,
        FUN = median)
    data[["roundedMedian"]] <- round(data[[medianCol]], digits = digits)
    geom_label(
        data = data,
        mapping = aes_string(label = "roundedMedian"),
        alpha = qcLabelAlpha,
        color = qcLabelColor,
        fill = qcLabelFill,
        fontface = qcLabelFontface,
        label.padding = qcLabelPadding,
        label.size = qcLabelSize,
        show.legend = FALSE)
}
