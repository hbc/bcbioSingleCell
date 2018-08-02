#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @name plotUMIsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @inheritParams general
#' @param point `string`. Label either the "`knee`" or "`inflection`" points per
#'   sample. To disable, use "`none`". Requires `geom = "ecdf"`.
#'
#' @examples
#' plotUMIsPerCell(indrops_small, geom = "violin")
#' plotUMIsPerCell(indrops_small, geom = "ridgeline")
#' plotUMIsPerCell(indrops_small, geom = "ecdf")
#' plotUMIsPerCell(indrops_small, geom = "histogram")
#' plotUMIsPerCell(indrops_small, geom = "boxplot")
NULL



# Methods ======================================================================
#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        interestingGroups,
        min = 0L,
        max = Inf,
        point = c("none", "inflection", "knee"),
        trans = "log10",
        color = getOption("bcbio.discrete.color", NULL),
        fill = getOption("bcbio.discrete.fill", NULL),
        title = "UMIs per cell"
    ) {
        geom <- match.arg(geom)
        point <- match.arg(point)
        assertIsAStringOrNULL(title)

        # Override interestingGroups if labeling points
        if (point != "none") {
            interestingGroups <- "sampleName"
        }

        p <- .plotQCMetric(
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

        # Calculate barcode ranks and label inflection or knee points
        if (point != "none") {
            # Require ecdf geom for now
            assert_are_identical(geom, "ecdf")

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

            assert_are_identical(
                x = names(sampleNames),
                y = names(points)
            )

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
                    bcbio_geom_label_repel(
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
)
