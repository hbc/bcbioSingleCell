#' Plot Barcode Ranks
#'
#' @name plotBarcodeRanks
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Passthrough arguments to [DropletUtils::barcodeRanks()].
#'
#' @seealso [DropletUtils::barcodeRanks()].
#'
#' @examples
#' # SingleCellExperiment ====
#' plotBarcodeRanks(cellranger_small)
#'
#' # seurat ====
#' plotBarcodeRanks(seurat_small)
NULL



# Methods ======================================================================
#' @rdname plotBarcodeRanks
#' @export
setMethod(
    "plotBarcodeRanks",
    signature("SingleCellExperiment"),
    function(object, ...) {
        ranksPerSample <- barcodeRanksPerSample(object, ...)
        sampleNames <- sampleData(object) %>%
            .[names(ranksPerSample), "sampleName", drop = TRUE] %>%
            as.character()

        plotlist <- mapply(
            sampleName = sampleNames,
            ranks = ranksPerSample,
            FUN = function(sampleName, ranks) {
                data <- cbind(
                    "rank" = ranks[["rank"]],
                    "total" = ranks[["total"]],  # nUMI
                    "fitted" = ranks[["fitted"]]
                ) %>%
                    as("tbl_df")

                p <- ggplot(data = data) +
                    geom_point(
                        mapping = aes_string(
                            x = "rank",
                            y = "total"
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
                    mapping = aes_string(
                        x = "rank",
                        y = "fitted"
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
                labelData <- rbind(
                    data[
                        which.min(abs(data[["total"]] - ranks[["knee"]]))
                        ,
                    ],
                    data[
                        which.min(abs(data[["total"]] - ranks[["inflection"]]))
                        ,
                    ]
                )
                labelData[["label"]] <- c(
                    paste("knee", "=", ranks[["knee"]]),
                    paste("inflection", "=", ranks[["inflection"]])
                )
                p +
                    bcbio_geom_label_repel(
                        data = labelData,
                        mapping = aes_string(
                            x = "rank",
                            y = "total",
                            label = "label"
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
)



#' @rdname plotBarcodeRanks
#' @export
setMethod(
    "plotBarcodeRanks",
    signature("seurat"),
    getMethod("plotBarcodeRanks", "SingleCellExperiment")
)
