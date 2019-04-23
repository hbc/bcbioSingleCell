#' @name plotBarcodeRanks
#' @author Michael Steinbaugh
#' @include barcodeRanksPerSample-methods.R
#' @inherit bioverbs::plotBarcodeRanks
#' @inherit barcodeRanksPerSample
#' @examples
#' data(indrops)
#' plotBarcodeRanks(indrops)
NULL



#' @rdname plotBarcodeRanks
#' @name plotBarcodeRanks
#' @importFrom bioverbs plotBarcodeRanks
#' @export
NULL



plotBarcodeRanks.bcbioSingleCell <-  # nolint
    function(
        object,
        colors = set_names(
            x = synesthesia(n = 3L),
            value = c("knee", "inflection", "fitline")
        )
    ) {
        ranksPerSample <- do.call(
            what = barcodeRanksPerSample,
            args = matchArgsToDoCall(args = list(object = object))
        )

        sampleData <- sampleData(object)
        if (is.null(sampleData)) {
            sampleNames <- "unknown"
        } else {
            sampleNames <- sampleData(object) %>%
                .[names(ranksPerSample), "sampleName", drop = TRUE] %>%
                as.character()
        }

        plotlist <- mapply(
            sampleName = sampleNames,
            ranks = ranksPerSample,
            FUN = function(sampleName, ranks) {
                data <- tibble(
                    rank = ranks[["rank"]],
                    total = ranks[["total"]],  # nUMI
                    fitted = ranks[["fitted"]]
                )

                p <- ggplot(data = data) +
                    geom_point(
                        mapping = aes(
                            x = !!sym("rank"),
                            y = !!sym("total")
                        )
                    ) +
                    scale_x_continuous(trans = "log10") +
                    scale_y_continuous(trans = "log10") +
                    labs(
                        title = sampleName,
                        y = "UMIs per cell"
                    )

                # Include the fit line (smooth.spline)
                p <- p + geom_line(
                    data = filter(data, !is.na(!!sym("fitted"))),
                    mapping = aes(
                        x = !!sym("rank"),
                        y = !!sym("fitted")
                    ),
                    colour = colors[["fitline"]],
                    size = 1L
                )

                p <- p +
                    geom_hline(
                        colour = colors[["knee"]],
                        linetype = "dashed",
                        yintercept = ranks[["knee"]]
                    ) +
                    geom_hline(
                        colour = colors[["inflection"]],
                        linetype = "dashed",
                        yintercept = ranks[["inflection"]]
                    )

                # Label the knee and inflection points more clearly
                knee <- which.min(abs(
                    data[["total"]] - ranks[["knee"]]
                ))
                inflection <- which.min(abs(
                    data[["total"]] - ranks[["inflection"]]
                ))
                labelData <- data[c(knee, inflection), ]
                labelData[["label"]] <- c(
                    paste("knee", "=", ranks[["knee"]]),
                    paste("inflection", "=", ranks[["inflection"]])
                )
                p +
                    acid_geom_label_repel(
                        data = labelData,
                        mapping = aes(
                            x = !!sym("rank"),
                            y = !!sym("total"),
                            label = !!sym("label")
                        ),
                        colour = colors[c("knee", "inflection")]
                    )
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )

        # Sort the plots by sample name
        plotlist <- plotlist[sort(names(plotlist))]

        plot_grid(plotlist = plotlist)
    }

formals(plotBarcodeRanks.bcbioSingleCell) <-
    formals(barcodeRanksPerSample.bcbioSingleCell)



#' @rdname plotBarcodeRanks
#' @export
setMethod(
    f = "plotBarcodeRanks",
    signature = signature("bcbioSingleCell"),
    definition = plotBarcodeRanks.bcbioSingleCell
)
