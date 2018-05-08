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
#' @param point Label either the "`knee`" or "`inflection`" points per sample.
#'   To disable, use "`none`". Requires `geom = "ecdf"`.
#' @param label Label the points per sample on the plot. Only applies when
#'   `point` argument is defined.
#'
#' @examples
#' # SingleCellExperiment ====
#' # Label the inflection point per sample
#' plotUMIsPerCell(cellranger_small, geom = "ecdf", point = "inflection")
#'
#' # Label the knee point per sample
#' plotUMIsPerCell(cellranger_small, geom = "ecdf", point = "knee")
#'
#' # Alternate geoms
#' plotUMIsPerCell(cellranger_small, geom = "histogram")
#' plotUMIsPerCell(cellranger_small, geom = "ridgeline")
#' plotUMIsPerCell(cellranger_small, geom = "violin")
#' plotUMIsPerCell(cellranger_small, geom = "boxplot")
#'
#' # seurat ====
#' plotUMIsPerCell(seurat_small, geom = "ecdf", point = "inflection")
NULL



# Methods ======================================================================
#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("histogram", "ecdf", "ridgeline", "violin", "boxplot"),
        interestingGroups,
        min = 0L,
        point = c("none", "inflection", "knee"),
        label = TRUE,
        trans = "log10",
        color = scale_color_hue(),
        fill = scale_fill_hue(),
        title = "UMIs per cell"
    ) {
        geom <- match.arg(geom)
        point <- match.arg(point)
        assert_is_a_bool(label)
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

            ranks <- barcodeRanksPerSample(object)
            # Inflection or knee points per sample
            points <- lapply(seq_along(ranks), function(x) {
                ranks[[x]][[point]]
            })
            points <- unlist(points)
            names(points) <- names(ranks)

            sampleNames <- sampleNames(object)

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
                    "x" = points,
                    "y" = freq,
                    "label" = paste0(sampleNames, " (", points, ")"),
                    "sampleName" = sampleNames
                )
                p <- p +
                    geom_point(
                        data = pointData,
                        mapping = aes_string(
                            x = "x",
                            y = "y",
                            color = "sampleName"
                        ),
                        size = 5L,
                        show.legend = FALSE
                    )
                if (isTRUE(label)) {
                    p <- p +
                        bcbio_geom_label_repel(
                            data = pointData,
                            mapping = aes_string(
                                x = "x",
                                y = "y",
                                label = "label",
                                color = "sampleName"
                            )
                        )
                }
            }
        }

        p <- p + labs(title = title)

        p
    }
)



#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("seurat"),
    getMethod("plotUMIsPerCell", "SingleCellExperiment")
)
