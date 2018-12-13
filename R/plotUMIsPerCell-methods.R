#' @name plotUMIsPerCell
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#' @inherit basejump::plotUMIsPerCell
#' @inheritParams basejump::params

#' @param point `character(1)`. Label either the "`knee`" or "`inflection`" points per
#'   sample. To disable, use "`none`". Requires `geom = "ecdf"`.
#'
#' @examples
#' data(indrops)
#' plotUMIsPerCell(indrops, geom = "violin")
#' plotUMIsPerCell(indrops, geom = "ridgeline")
#' plotUMIsPerCell(indrops, geom = "ecdf")
#' plotUMIsPerCell(indrops, geom = "histogram")
#' plotUMIsPerCell(indrops, geom = "boxplot")
NULL



#' @importFrom basejump plotUMIsPerCell
#' @aliases NULL
#' @export
basejump::plotUMIsPerCell



plotUMIsPerCell.SingleCellExperiment <-  # nolint
    function(
        object,
        geom,
        interestingGroups = NULL,
        min = 0L,
        max = Inf,
        point = c("none", "inflection", "knee"),
        trans = "log10",
        color,
        fill,
        title = "UMIs per cell"
    ) {
        assert(isString(title) || is.null(title))
        geom <- match.arg(geom)
        point <- match.arg(point)

        # Override interestingGroups if labeling points
        if (point != "none") {
            interestingGroups <- "sampleName"
        }

        p <- do.call(
            what = .plotQCMetric,
            args = list(
                object = object,
                metricCol = "nUMI",
                geom = geom,
                interestingGroups = interestingGroups,
                min = min,
                max = max,
                trans = trans,
                color = color,
                fill = fill
            )
        )

        # Calculate barcode ranks and label inflection or knee points.
        if (point != "none") {
            # Require ecdf geom for now
            assert(identical(geom, "ecdf"))

            if (length(title)) {
                p <- p + labs(subtitle = paste(point, "point per sample"))
            }

            sampleNames <- sampleNames(object)

            ranks <- barcodeRanksPerSample(object)
            # Inflection or knee points per sample
            points <- lapply(seq_along(ranks), function(x) {
                ranks[[x]][[point]]
            })
            points <- unlist(points)
            names(points) <- names(ranks)

            assert(identical(names(sampleNames), names(points)))

            if (geom == "ecdf") {
                # Calculate the y-intercept per sample
                freq <- mapply(
                    sampleID = names(points),
                    point = points,
                    MoreArgs = list(metrics = metrics(object)),
                    FUN = function(metrics, sampleID, point) {
                        nUMI <- metrics[
                            metrics[["sampleID"]] == sampleID,
                            "nUMI",
                            drop = TRUE
                        ]
                        e <- ecdf(sort(nUMI))
                        e(point)
                    },
                    SIMPLIFY = TRUE,
                    USE.NAMES = TRUE
                )
                pointData <- data.frame(
                    x = points,
                    y = freq,
                    label = paste0(sampleNames, " (", points, ")"),
                    sampleName = sampleNames
                )
                p <- p +
                    geom_point(
                        data = pointData,
                        mapping = aes(
                            x = !!sym("x"),
                            y = !!sym("y"),
                            color = !!sym("sampleName")
                        ),
                        size = 5L,
                        show.legend = FALSE
                    ) +
                    basejump_geom_label_repel(
                        data = pointData,
                        mapping = aes(
                            x = !!sym("x"),
                            y = !!sym("y"),
                            label = !!sym("label"),
                            color = !!sym("sampleName")
                        )
                    )
            }
        }

        p <- p + labs(title = title)

        p
    }

formals(plotUMIsPerCell.SingleCellExperiment)[["color"]] <-
    formalsList[["color.discrete"]]
formals(plotUMIsPerCell.SingleCellExperiment)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotUMIsPerCell.SingleCellExperiment)[["geom"]] <- geom



#' @rdname plotUMIsPerCell
#' @export
setMethod(
    f = "plotUMIsPerCell",
    signature = signature("SingleCellExperiment"),
    definition = plotUMIsPerCell.SingleCellExperiment
)
