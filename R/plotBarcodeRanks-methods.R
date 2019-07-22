#' @name plotBarcodeRanks
#' @author Michael Steinbaugh
#' @include barcodeRanksPerSample-methods.R
#' @inherit bioverbs::plotBarcodeRanks
#' @inherit barcodeRanksPerSample
#'
#' @param colors `character(3)`.
#'   Character vector denoting `fitline`, `inflection`, and `knee` point colors.
#'   Must pass in color names or hexadecimal values.
#' @param ... Additional arguments.
#'
#' @examples
#' data(indrops)
#' if (packageVersion("DropletUtils") >= "1.4") {
#'     plotBarcodeRanks(indrops)
#' }
NULL



#' @rdname plotBarcodeRanks
#' @name plotBarcodeRanks
#' @importFrom bioverbs plotBarcodeRanks
#' @usage plotBarcodeRanks(object, ...)
#' @export
NULL



plotBarcodeRanks.bcbioSingleCell <-  # nolint
    function(
        object,
        colors = c(
            fitline = "dodgerblue",
            inflection = "darkorchid3",
            knee = "darkorange2"
        )
    ) {
        validObject(object)
        assert(
            isCharacter(colors),
            areSetEqual(
                x = names(colors),
                y = c("fitline", "inflection", "knee")
            )
        )

        ranksPerSample <- do.call(
            what = barcodeRanksPerSample,
            args = matchArgsToDoCall(
                args = list(object = object),
                removeFormals = "colors"
            )
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
                data <- as_tibble(ranks, rownames = NULL)
                inflection <- metadata(ranks)[["inflection"]]
                knee <- metadata(ranks)[["knee"]]

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

                ## Include the fit line (smooth.spline)
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
                        yintercept = knee
                    ) +
                    geom_hline(
                        colour = colors[["inflection"]],
                        linetype = "dashed",
                        yintercept = inflection
                    )

                ## Label the knee and inflection points more clearly
                knee <- which.min(abs(data[["total"]] - knee))
                inflection <- which.min(abs(data[["total"]] - inflection))
                labelData <- data[c(knee, inflection), , drop = FALSE]
                labelData[["label"]] <- c(
                    paste("knee", "=", knee),
                    paste("inflection", "=", inflection)
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

        ## Sort the plots by sample name
        plotlist <- plotlist[sort(names(plotlist))]

        plot_grid(plotlist = plotlist)
    }

f1 <- formals(plotBarcodeRanks.bcbioSingleCell)
f2 <- formals(barcodeRanksPerSample.bcbioSingleCell)
f2 <- f2[setdiff(names(f2), names(f1))]
f <- c(f1, f2)
formals(plotBarcodeRanks.bcbioSingleCell) <- f



#' @rdname plotBarcodeRanks
#' @export
setMethod(
    f = "plotBarcodeRanks",
    signature = signature("bcbioSingleCell"),
    definition = plotBarcodeRanks.bcbioSingleCell
)
