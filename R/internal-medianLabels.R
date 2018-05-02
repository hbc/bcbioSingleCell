#' Add Median Labels to a Plot
#'
#' For example, `medianCol` can be `nGene`. Median values are always calculated
#' by sample (`sampleName`).
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @return ggplot2 `geom_label`.
.medianLabels <- function(data, medianCol, digits = 0L) {
    assert_is_data.frame(data)
    assert_is_a_string(medianCol)
    assert_is_subset(medianCol, colnames(data))
    assert_is_an_integer(digits)

    data <- aggregate(
        formula = as.formula(paste(medianCol, "sampleName", sep = " ~ ")),
        data = data,
        FUN = median
    )
    data[["roundedMedian"]] <- round(data[[medianCol]], digits = digits)

    # Add `aggregate` column for facet wrapping, if necessary
    if ("aggregate" %in% colnames(data)) {
        sampleFacet <- data[, c("sampleName", "aggregate")] %>%
        unique()
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
