#' General Arguments
#'
#' @name general
#' @keywords internal
#'
#' @param color `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot2 color
#'   scale. Must supply discrete values. When set to `NULL`, the default ggplot2
#'   color palette will be used. If manual color definitions are desired, we
#'   recommend using [ggplot2::scale_color_manual()].
#'   To set the discrete color palette globally, use
#'   `options(bcbio.discrete.color = scale_color_viridis_d())`.
#' @param fill `ggproto`/`ScaleDiscrete` or `NULL`. Desired ggplot2 fill scale.
#'   Must supply discrete values. When set to `NULL`, the default ggplot2 color
#'   palette will be used. If manual color definitions are desired, we recommend
#'   using [ggplot2::scale_fill_manual()].
#'   To set the discrete fill palette globally, use
#'   `options(bcbio.discrete.fill = scale_fill_viridis_d())`.
#' @param geom `string`. Plot type. Uses [match.arg()] and defaults to the first
#'   argument in the `character` vector.
#' @param interestingGroups `character` or `NULL`. Character vector of
#'   interesting groups. Must be formatted in camel case and intersect with
#'   [sampleData()] colnames.
#' @param label `boolean`. Overlay a cluster identitiy label on the plot.
#' @param labelSize `scalar integer`. Size of the text label.
#' @param legend `boolean`. Include plot legend.
#' @param min `scalar numeric`. Recommended minimum value cutoff.
#' @param max `scalar numeric`. Recommended maximum value cutoff.
#' @param object Object.
#' @param prefilter `boolean`. Apply prefiltering to remove zero count genes.
#' @param pipeline `string`. Pipeline used to generate the samples.
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
