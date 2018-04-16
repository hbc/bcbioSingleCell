.labelBarcodeRanks <- function(
    p,
    object,
    geom,
    point = c("knee", "inflection")
) {
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
    sampleNames <- as.character(sampleData(object)[, "sampleName"])

    points <- lapply(seq_along(ranks), function(x) {
        ranks[[x]][[point]]
    })
    points <- unlist(points)
    names(points) <- names(ranks)

    labels <- paste(sampleNames, points, sep = " = ")

    p <- p + labs(
        subtitle = paste(labels, collapse = "\n")
    )

    p
}
