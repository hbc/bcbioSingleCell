.labelBarcodeRanks <- function(
    p,
    object,
    geom,
    point = c("knee", "inflection")
) {
    .Deprecated(".labelBarcodeRanksPerSample")
    point <- match.arg(point)
    ranks <- barcodeRanks(object)
    inflection <- ranks[["inflection"]]
    knee <- ranks[["knee"]]

    if (point == "inflection") {
        if (geom %in% c("boxplot", "violin")) {
            xintercept <- NULL
            yintercept <- inflection
        } else {
            xintercept <- inflection
            yintercept <- NULL
        }
        p <- p +
            .qcCutoffLine(
                xintercept = xintercept,
                yintercept = yintercept,
                color = inflectionColor
            )
        subtitle <- paste("inflection", inflection, sep = " = ")
    }

    if (point == "knee") {
        if (geom %in% c("boxplot", "violin")) {
            xintercept <- NULL
            yintercept <- knee
        } else {
            xintercept <- knee
            yintercept <- NULL
        }
        p <- p +
            .qcCutoffLine(
                xintercept = xintercept,
                yintercept = yintercept,
                color = kneeColor
            )
        subtitle <- paste("knee", knee, sep = " = ")
    }

    p <- p + labs(subtitle = subtitle)

    p
}



.labelBarcodeRanksPerSample <- function(
    p,
    object,
    geom,
    point = c("knee", "inflection")
) {
    point <- match.arg(point)
    ranks <- barcodeRanksPerSample(object)

    points <- lapply(seq_along(ranks), function(x) {
        ranks[[x]][[point]]
    })
    points <- unlist(points)
    names(points) <- names(ranks)

    # TODO Use `ecdf()` to calculate y point instead
    # This doesn't seem to be fitting quite right, work on the barcodeRanks method
    sampleName <- sampleData(object) %>%
        .[names(points), "sampleName"] %>%
        as.character()
    labelData <- data.frame(
        "x" = points,
        "y" = 1L,
        "label" = paste0(sampleName, " (", points, ")"),
        "sampleName" = sampleName
    )

    # FIXME Map colors to sample names
    p <- p +
        geom_label_repel(
            data = labelData,
            mapping = aes_string(
                x = "x",
                y = "y",
                label = "label",
                color = "sampleName"
            ),
            fill = "white",
            show.legend = FALSE
        )

    p
}
