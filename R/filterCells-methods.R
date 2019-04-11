#' @name filterCells
#' @author Michael Steinbaugh
#' @inherit bioverbs::filterCells
#' @inheritParams basejump::params
#'
#' @details
#' Apply gene detection, novelty score, and mitochondrial abundance cutoffs to
#' cellular barcodes. By default we recommend applying the same filtering cutoff
#' to all samples. The filtering parameters now support per-sample cutoffs,
#' defined using a named `numeric` vector. When matching per sample, be sure to
#' use the `sampleNames` return values (i.e. the `sampleName` column in
#' `sampleData`).
#'
#' @param nCells `integer(1)`.
#'   Expected number of cells per sample.
#' @param minUMIs `integer(1)`.
#'   Minimum number of UMI disambiguated counts per cell.
#' @param maxUMIs `integer(1)`.
#'   Maximum number of UMI disambiguated counts per cell.
#' @param minGenes `integer(1)`.
#'   Minimum number of genes detected.
#' @param maxGenes `integer(1)`.
#'   Maximum number of genes detected.
#' @param minNovelty `integer(1)` (`0`-`1`).
#'   Minimum novelty score (log10 genes per UMI).
#' @param maxMitoRatio `integer(1)` (`0`-`1`).
#'   Maximum relative mitochondrial abundance.
#' @param minCellsPerGene `integer(1)`.
#'   Include genes with non-zero expression in at least this many cells.
#'
#' @return `bcbioSingleCell`.
#' Filtering information gets slotted into [`metadata()`][S4Vectors::metadata]
#' as `filterCells` and `filterParams`.
#'
#' @examples
#' data(indrops)
#'
#' object <- indrops
#' show(object)
#'
#' x <- filterCells(object)
#' show(x)
#' metadata(x)$filterParams
#'
#' ## Per sample cutoffs.
#' sampleNames(object)
#' x <- filterCells(
#'     object = object,
#'     minUMIs = c(rep_1 = 100)
#' )
#' show(x)
#' metadata(x)$filterParams
NULL



#' @rdname filterCells
#' @name filterCells
#' @importFrom bioverbs filterCells
#' @export
NULL



.isFiltered <- function(object) {
    if (!is.null(metadata(object)[["filterParams"]])) {
        TRUE
    } else {
        FALSE
    }
}



.paddedCount <- function(x, width = 8L) {
    str_pad(x, width = width, pad = " ")
}



filterCells.bcbioSingleCell <-  # nolint
    function(
        object,
        nCells = Inf,
        minUMIs = 1L,
        maxUMIs = Inf,
        minGenes = 1L,
        maxGenes = Inf,
        minNovelty = 0L,
        maxMitoRatio = 1L,
        minCellsPerGene = 1L
    ) {
        validObject(object)

        originalDim <- dim(object)
        sampleNames <- sampleNames(object)
        # Coercing to tibble here for dplyr operations.
        colData <- as(object = colData(object), Class = "tbl_df")

        # Parameter integrity checks -------------------------------------------
        assert(
            # Expected nCells per sample.
            is.numeric(nCells),
            all(isPositive(nCells)),

            # minUMIs (see below)
            isAny(minUMIs, c("numeric", "character")),

            # maxUMIs
            is.numeric(maxUMIs),
            all(isPositive(maxUMIs)),

            # minGenes
            is.numeric(minGenes),
            all(isPositive(minGenes)),

            # maxGenes
            is.numeric(maxGenes),
            all(isNonNegative(maxGenes)),

            # minNovelty
            is.numeric(minNovelty),
            all(isInRange(minNovelty, lower = 0L, upper = 1L)),

            # maxMitoRatio
            is.numeric(maxMitoRatio),
            # Don't allow the user to set at 0.
            all(isInLeftOpenRange(maxMitoRatio, lower = 0L, upper = 1L)),

            # minCellsPerGene
            is.numeric(minCellsPerGene),
            # Don't allow genes with all zero counts.
            all(isPositive(minCellsPerGene))
        )

        # minUMIs supports barcode ranks filtering.
        if (is.character(minUMIs)) {
            assert(
                isString(minUMIs),
                isSubset(minUMIs, c("inflection", "knee"))
            )
        } else {
            assert(
                is.numeric(minUMIs),
                all(isPositive(minUMIs))
            )
        }

        params <- list(
            nCells = nCells,
            minUMIs = minUMIs,
            maxUMIs = maxUMIs,
            minGenes = minGenes,
            maxGenes = maxGenes,
            minNovelty = minNovelty,
            maxMitoRatio = maxMitoRatio,
            minCellsPerGene = minCellsPerGene
        )

        # Filter low quality cells ---------------------------------------------
        summaryCells <- character()
        summaryCells[["prefilter"]] <- paste(
            paste(.paddedCount(ncol(object)), "cells"),
            "prefilter",
            sep = " | "
        )

        # minUMIs
        if (isString(minUMIs)) {
            ranks <- barcodeRanksPerSample(object)
            minUMIs <- vapply(
                X = ranks,
                FUN = function(x) {
                    as.integer(x[[minUMIs]])
                },
                FUN.VALUE = integer(1L)
            )
            names(minUMIs) <- sampleNames
            minUMIs <- minUMIs[sort(names(minUMIs))]
        }
        if (!is.null(names(minUMIs))) {
            assert(areSetEqual(names(minUMIs), sampleNames))
            message(paste(
                "minUMIs: per sample mode",
                printString(minUMIs),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(minUMIs),
                cutoff = minUMIs,
                FUN = function(sample, cutoff) {
                    filter(
                        colData,
                        !!sym("sampleName") == !!sample,
                        !!sym("nUMI") >= !!cutoff
                    )
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            colData <- bind_rows(list)
        } else {
            colData <- filter(colData, !!sym("nUMI") >= !!minUMIs)
        }
        if (!nrow(colData)) {
            stop("No cells passed `minUMIs` cutoff")
        }
        summaryCells[["minUMIs"]] <- paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            paste("minUMIs", ">=", min(minUMIs)),
            sep = " | "
        )

        # maxUMIs
        if (!is.null(names(maxUMIs))) {
            assert(areSetEqual(names(maxUMIs), sampleNames))
            message(paste(
                "maxUMIs: per sample mode",
                printString(maxUMIs),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(maxUMIs),
                cutoff = maxUMIs,
                FUN = function(sample, cutoff) {
                    filter(
                        colData,
                        !!sym("sampleName") == !!sample,
                        !!sym("nUMI") <= !!cutoff
                    )
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            colData <- bind_rows(list)
        } else {
            colData <- filter(colData, !!sym("nUMI") <= !!maxUMIs)
        }
        if (!nrow(colData)) {
            stop("No cells passed `maxUMIs` cutoff")
        }
        summaryCells[["maxUMIs"]] <- paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            paste("maxUMIs", "<=", max(maxUMIs)),
            sep = " | "
        )

        # minGenes
        if (!is.null(names(minGenes))) {
            assert(areSetEqual(names(minGenes), sampleNames))
            message(paste(
                "minGenes: per sample mode",
                printString(minGenes),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(minGenes),
                cutoff = minGenes,
                FUN = function(sample, cutoff) {
                    filter(
                        colData,
                        !!sym("sampleName") == !!sample,
                        !!sym("nGene") >= !!cutoff
                    )
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            colData <- bind_rows(list)
        } else {
            colData <- filter(colData, !!sym("nGene") >= !!minGenes)
        }
        if (!nrow(colData)) {
            stop("No cells passed `minGenes` cutoff")
        }
        summaryCells[["minGenes"]] <- paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            paste("minGenes", ">=", min(minGenes)),
            sep = " | "
        )

        # maxGenes
        if (!is.null(names(maxGenes))) {
            assert(areSetEqual(names(maxGenes), sampleNames))
            message(paste(
                "maxGenes: per sample mode",
                printString(maxGenes),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(maxGenes),
                cutoff = maxGenes,
                FUN = function(sample, cutoff) {
                    filter(
                        colData,
                        !!sym("sampleName") == !!sample,
                        !!sym("nGene") <= !!cutoff
                    )
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            colData <- bind_rows(list)
        } else {
            colData <- filter(colData, !!sym("nGene") <= !!maxGenes)
        }
        if (!nrow(colData)) {
            stop("No cells passed `maxGenes` cutoff")
        }
        summaryCells[["maxGenes"]] <- paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            paste("maxGenes", "<=", max(maxGenes)),
            sep = " | "
        )

        # minNovelty
        if (!is.null(names(minNovelty))) {
            assert(areSetEqual(names(minNovelty), sampleNames))
            message(paste(
                "minNovelty: per sample mode",
                printString(minNovelty),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(minNovelty),
                cutoff = minNovelty,
                FUN = function(sample, cutoff) {
                    filter(
                        colData,
                        !!sym("sampleName") == !!sample,
                        !!sym("log10GenesPerUMI") >= !!cutoff)
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            colData <- bind_rows(list)
        } else {
            colData <- filter(
                colData,
                !!sym("log10GenesPerUMI") >= !!minNovelty
            )
        }
        if (!nrow(colData)) {
            stop("No cells passed `minNovelty` cutoff")
        }
        summaryCells[["minNovelty"]] <- paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            paste("minNovelty", "<=", min(minNovelty)),
            sep = " | "
        )

        # maxMitoRatio
        if (!is.null(names(maxMitoRatio))) {
            assert(areSetEqual(names(maxMitoRatio), sampleNames))
            message(paste(
                "maxMitoRatio: per sample mode",
                printString(maxMitoRatio),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(maxMitoRatio),
                cutoff = maxMitoRatio,
                FUN = function(sample, cutoff) {
                    filter(
                        colData,
                        !!sym("sampleName") == !!sample,
                        !!sym("mitoRatio") <= !!cutoff
                    )
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            colData <- bind_rows(list)
        } else {
            colData <- filter(colData, !!sym("mitoRatio") <= !!maxMitoRatio)
        }
        if (!nrow(colData)) {
            stop("No cells passed `maxMitoRatio` cutoff")
        }
        summaryCells[["maxMitoRatio"]] <- paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            paste("maxMitoRatio", "<=", max(maxMitoRatio)),
            sep = " | "
        )

        # Expected nCells per sample (filtered by top nUMI)
        if (nCells < Inf) {
            colData <- colData %>%
                group_by(!!sym("sampleID")) %>%
                arrange(desc(!!sym("nUMI")), .by_group = TRUE) %>%
                slice(seq_len(nCells))
        }

        if (!nrow(colData)) {
            stop("No cells passed `nCells` cutoff")
        }
        summaryCells[["nCells"]] <- paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            paste("nCells", "==", nCells),
            sep = " | "
        )

        # Now coerce back to DataFrame from tibble.
        colData <- as(colData, "DataFrame")
        assert(hasRownames(colData))
        cells <- sort(rownames(colData))
        assert(isSubset(cells, colnames(object)))
        object <- object[, cells, drop = FALSE]

        # Filter low quality genes ---------------------------------------------
        summaryGenes <- character()
        summaryGenes[["prefilter"]] <- paste(
            paste(.paddedCount(nrow(object)), "genes"),
            "prefilter",
            sep = " | "
        )
        if (minCellsPerGene > 0L) {
            nonzero <- counts(object) > 0L
            keep <- rowSums(nonzero) >= minCellsPerGene
            genes <- names(keep)[keep]
        } else {
            genes <- rownames(object)
        }
        if (!length(genes)) {
            stop("No genes passed `minCellsPerGene` cutoff")
        }
        summaryGenes[["minCellsPerGene"]] <- paste(
            paste(.paddedCount(length(genes)), "genes"),
            paste("minCellsPerGene", ">=", as.character(minCellsPerGene)),
            sep = " | "
        )
        genes <- sort(genes)
        assert(isSubset(genes, rownames(object)))
        object <- object[genes, , drop = FALSE]

        # Summary --------------------------------------------------------------
        summaryParams <- paste("  -", c(
            paste(">=", min(minUMIs), "UMIs per cell"),
            paste("<=", max(maxUMIs), "UMIs per cell"),
            paste(">=", min(minGenes), "genes per cell"),
            paste("<=", max(maxGenes), "genes per cell"),
            paste(">=", min(minNovelty), "novelty score"),
            paste("<=", max(maxMitoRatio), "mitochondrial abundance"),
            paste("==", nCells, "cells per sample"),
            paste(">=", min(minCellsPerGene), "cells per gene")
        ))
        summary <- c(
            "Parameters:",
            summaryParams,
            separatorBar,
            "Cells:",
            as.character(summaryCells),
            paste(
                .paddedCount(dim(object)[[2L]]), "of", originalDim[[2L]],
                "cells passed filtering",
                paste0(
                    "(", percent(dim(object)[[2L]] / originalDim[[2L]]), ")"
                )
            ),
            # Number of cells per sample.
            printString(table(colData[["sampleName"]])),
            separatorBar,
            "Genes:",
            as.character(summaryGenes),
            paste(
                .paddedCount(dim(object)[[1L]]), "of", originalDim[[1L]],
                "genes passed filtering",
                paste0(
                    "(", percent(dim(object)[[1L]] / originalDim[[1L]]), ")"
                )
            )
        )
        message(paste(summary, collapse = "\n"))

        # Metadata -------------------------------------------------------------
        metadata(object)[["cellularBarcodes"]] <- NULL
        metadata(object)[["filterCells"]] <- cells
        metadata(object)[["filterGenes"]] <- genes
        metadata(object)[["filterParams"]] <- params
        metadata(object)[["filterSummary"]] <- summary

        object
    }



#' @rdname filterCells
#' @export
setMethod(
    f = "filterCells",
    signature = signature("bcbioSingleCell"),
    definition = filterCells.bcbioSingleCell
)
