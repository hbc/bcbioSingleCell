# FIXME Use the formals for plotBarcodeRanks



#' Plot Barcode Ranks
#'
#'
#' @name plotBarcodeRanks
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inherit barcodeRanksPerSample
#'
#' @return `ggplot` grid.
#'
#' @examples
#' plotBarcodeRanks(indrops_small)
NULL



.plotBarcodeRanks <- function(object) {
    ranksPerSample <- do.call(
        what = barcodeRanksPerSample,
        args = setArgsToDoCall(
            args = list(object = object),
            call = matchCall()
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
            data <- cbind(
                rank = ranks[["rank"]],
                total = ranks[["total"]],  # nUMI
                fitted = ranks[["fitted"]]
            ) %>%
                as("tbl_df")

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
                color = "red",
                size = 1L
            )

            # Knee and inflection points
            colors <- c("dodgerblue", "forestgreen")
            p <- p +
                geom_hline(
                    color = colors[[1L]],
                    linetype = "dashed",
                    yintercept = ranks[["knee"]]
                ) +
                geom_hline(
                    color = colors[[2L]],
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
                basejump_geom_label_repel(
                    data = labelData,
                    mapping = aes(
                        x = !!sym("rank"),
                        y = !!sym("total"),
                        label = !!sym("label")
                    ),
                    color = colors
                )
        },
        SIMPLIFY = FALSE,
        USE.NAMES = TRUE
    )

    # Sort the plots by sample name
    plotlist <- plotlist[sort(names(plotlist))]

    plot_grid(plotlist = plotlist)
}

formals(.plotBarcodeRanks) <- methodFormals(
    f = "barcodeRanksPerSample",
    signature = "SingleCellExperiment"
)



#' @rdname plotBarcodeRanks
#' @export
setMethod(
    f = "plotBarcodeRanks",
    signature = signature("SingleCellExperiment"),
    definition = .plotBarcodeRanks
)
