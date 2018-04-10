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
    color <- viridis(length(ranks))
    for (i in seq_along(ranks)) {
        if (point == "inflection") {
            inflection <- ranks[[i]][["inflection"]]
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
                    linetype = "solid",
                    alpha = 0.5,
                    size = 0.75,
                    color = color[[i]]
                )
        }
        if (point == "knee") {
            knee <- ranks[[i]][["knee"]]
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
                    linetype = "solid",
                    alpha = 0.5,
                    size = 0.75,
                    color = color[[i]]
                )
        }
    }

    p
}
