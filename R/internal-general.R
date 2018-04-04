#' General Arguments
#'
#' @name general
#' @keywords internal
#'
#' @param color Desired ggplot color scale. Must supply discrete values. When
#'   set to `NULL`, the default ggplot2 color palette will be used. If manual
#'   color definitions are desired, we recommend using
#'   [ggplot2::scale_color_manual()].
#' @param fill Desired ggplot fill scale. Must supply discrete values. When set
#'   to `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [ggplot2::scale_fill_manual()].
#' @param gene2symbol `data.frame` containing gene-to-symbol mappings. Columns
#'   must contain `geneID` and `geneName`.
#' @param dark Plot against a dark background using [basejump::midnightTheme()].
#' @param dir Output directory.
#' @param genes Gene identifiers. Must match the rownames of the object.
#' @param geom Plot type. Uses [match.arg()] internally and defaults to the
#'   first argument in the vector.
#' @param headerLevel R Markdown header level.
#' @param interestingGroups Character vector of interesting groups. Must be
#'   formatted in camel case and intersect with [sampleData()] colnames.
#' @param label Overlay a cluster identitiy label on the plot.
#' @param labelSize Size of the text label.
#' @param legend Include plot legend.
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#' @param object Object.
#' @param pipeline Pipeline used to generate the samples.
#' @param pointAlpha Alpha transparency level. Useful when there many cells in
#'   the dataset, and some cells can be masked.
#' @param pointsAsNumbers Plot the points as numbers (`TRUE`) or dots (`FALSE`).
#' @param pointSize Cell point size.
#' @param return Return type. Uses [base::match.arg()] internally and defaults
#'   to the first argument in the vector.
#' @param title Plot title.
#' @param trans Name of the axis scale transformation to apply. See
#'   `help("scale_x_continuous", "ggplot2")` for more information.
#' @param value Object to assign.
#' @param x Primary object.
#' @param y Secondary object.
#' @param ... Additional arguments.
#'
#' @return No value.
NULL
