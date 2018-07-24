#' General Arguments
#'
#' @name general
#' @keywords internal
#'
#' @param aspectRatio `scalar integer`. Aspect ratio.
#' @param color `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot color scale.
#'   Must supply discrete values. When set to `NULL`, the default ggplot2 color
#'   palette will be used. If manual color definitions are desired, we recommend
#'   using [ggplot2::scale_color_manual()].
#' @param fill `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot fill scale.
#'   Must supply discrete values. When set to `NULL`, the default ggplot2 color
#'   palette will be used. If manual color definitions are desired, we recommend
#'   using [ggplot2::scale_fill_manual()].
#' @param gene2symbol `data.frame`. Gene-to-symbol mappings. Columns must
#'   contain `geneID` and `geneName`.
#' @param dark `boolean`. Plot against a dark background using
#'   [basejump::theme_midnight()].
#' @param dir `string`. Output directory path.
#' @param expression `string`. Calculation to apply. Uses [match.arg()] and
#'   defaults to the first argument in the `character` vector.
#' @param genes `character`. Gene identifiers. Must match the rownames of the
#'   object.
#' @param geom `string`. Plot type. Uses [match.arg()] and defaults to the first
#'   argument in the `character` vector.
#' @param grid `boolean`. Show major grid lines but hide axis lines.
#' @param headerLevel `scalar integer` (`1`-`7`). R Markdown header level.
#' @param interestingGroups `character` or `NULL`. Character vector of
#'   interesting groups. Must be formatted in camel case and intersect with
#'   [sampleData()] colnames.
#' @param label `boolean`. Overlay a cluster identitiy label on the plot.
#' @param labelSize `scalar integer`. Size of the text label.
#' @param legend `boolean`. Include plot legend.
#' @param min `scalar numeric`. Recommended minimum value cutoff.
#' @param max `scalar numeric`. Recommended maximum value cutoff.
#' @param object Object.
#' @param pipeline `string`. Pipeline used to generate the samples.
#' @param pointAlpha `scalar numeric` (`0`-`1`). Alpha transparency level.
#'   Useful when there many cells in the dataset, and some cells can be masked.
#' @param pointsAsNumbers `boolean`. Plot the points as numbers (`TRUE`) or
#'   dots (`FALSE`).
#' @param pointSize `scalar numeric`. Cell point size.
#' @param reduction `string`. Dimensional reduction method to apply. Defaults to
#'   t-SNE ("`TSNE`") but UMAP is also supported ("`UMAP`").
#' @param return `string`. Return type. Uses [base::match.arg()] internally and
#'   defaults to the first argument in the `character` vector.
#' @param title `string` or `NULL`. Plot title.
#' @param trans `string`. Name of the axis scale transformation to apply. See
#'   `help("scale_x_continuous", "ggplot2")` for more information.
#' @param trendline `boolean`. Include trendline on plot.
#' @param value Object to assign.
#' @param x Primary object.
#' @param y Secondary object.
#' @param ... Additional arguments.
#'
#' @return No value.
NULL
