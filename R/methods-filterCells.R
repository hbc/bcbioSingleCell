#' Filter Cells
#'
#' Apply gene detection, novelty score, and mitochondrial abundance cutoffs to
#' cellular barcodes. By default we recommend applying the same filtering cutoff
#' to all samples. The filtering parameters now support per-sample cutoffs,
#' defined using a named `numeric` vector. When matching per sample, be sure to
#' use the [sampleNames()] return values (i.e. the `sampleName` column in
#' [sampleData()]).
#'
#' @name filterCells
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param minUMIs Minimum number of UMI disambiguated counts per cell.
#' @param maxUMIs Maximum number of UMI disambiguated counts per cell.
#' @param minGenes Minimum number of genes detected.
#' @param maxGenes Maximum number of genes detected.
#' @param minNovelty Minimum novelty score (log10 genes per UMI).
#' @param maxMitoRatio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param minCellsPerGene Include genes with non-zero expression in at least
#'   this many cells.
#'
#' @seealso [Seurat::CreateSeuratObject()].
#'
#' @return `bcbioSingleCell`, with filtering information slotted into
#'   [metadata()] as `filterCells` and `filterParams`.
#'
#' @examples
#' # bcbioSingleCell ====
#' show(indrops_small)
#' x <- filterCells(indrops_small, minGenes = 0L)
#' show(x)
#'
#' # SingleCellExperiment ====
#' show(cellranger_small)
#' x <- filterCells(cellranger_small, minNovelty = 0L)
#' show(x)
#'
#' # Per sample cutoffs
#' sampleNames(cellranger_small)
#' x <- filterCells(
#'     object = cellranger_small,
#'     minUMIs = c(
#'         distal = 100,
#'         proximal = 200
#'     )
#' )
#' x
NULL



# Constructors =================================================================
.paddedCount <- function(x, width = 8L) {
    str_pad(x, width = width, pad = " ")
}



# Methods ======================================================================
#' @rdname filterCells
#' @export
setMethod(
    "filterCells",
    signature("SingleCellExperiment"),
    function(
        object,
        minUMIs = 0L,
        maxUMIs = Inf,
        minGenes = 0L,
        maxGenes = Inf,
        minNovelty = 0L,
        maxMitoRatio = 1L,
        minCellsPerGene = 10L
    ) {
        validObject(object)
        sampleNames <- sampleNames(object)
        metrics <- metrics(object)

        # Parameter integrity checks ===========================================
        assert_is_any_of(minUMIs, c("numeric", "character"))
        if (is.character(minUMIs)) {
            assert_is_a_string(minUMIs)
            assert_is_subset(minUMIs, c("inflection", "knee"))
        }
        params <- list(
            minUMIs = minUMIs,
            maxUMIs = maxUMIs,
            minGenes = minGenes,
            maxGenes = maxGenes,
            minNovelty = minNovelty,
            maxMitoRatio = maxMitoRatio,
            minCellsPerGene = minCellsPerGene
        )

        # Filter low quality cells =============================================
        summaryCells <- character()
        summaryCells[["prefilter"]] <- paste(
            paste(.paddedCount(ncol(object)), "cells"),
            "prefilter",
            sep = " | "
        )

        # minUMIs --------------------------------------------------------------
        if (is_a_string(minUMIs)) {
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
            assert_are_set_equal(names(minUMIs), sampleNames)
            message(paste(
                "minUMIs: per sample mode",
                printString(minUMIs),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(minUMIs),
                cutoff = minUMIs,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleName"]] == sample, , drop = FALSE] %>%
                        .[.[["nUMI"]] >= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            metrics <- metrics %>%
                .[.[["nUMI"]] >= minUMIs, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `minUMIs` cutoff")
        }
        summaryCells[["minUMIs"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("minUMIs", ">=", min(minUMIs)),
            sep = " | "
        )

        # maxUMIs --------------------------------------------------------------
        if (!is.null(names(maxUMIs))) {
            assert_are_set_equal(names(maxUMIs), sampleNames)
            message(paste(
                "maxUMIs: per sample mode",
                printString(maxUMIs),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(maxUMIs),
                cutoff = maxUMIs,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleName"]] == sample, , drop = FALSE] %>%
                        .[.[["nUMI"]] <= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            metrics <- metrics %>%
                .[.[["nUMI"]] <= maxUMIs, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `maxUMIs` cutoff")
        }
        summaryCells[["maxUMIs"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("maxUMIs", "<=", max(maxUMIs)),
            sep = " | "
        )

        # minGenes -------------------------------------------------------------
        if (!is.null(names(minGenes))) {
            assert_are_set_equal(names(minGenes), sampleNames)
            message(paste(
                "minGenes: per sample mode",
                printString(minGenes),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(minGenes),
                cutoff = minGenes,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleName"]] == sample, , drop = FALSE] %>%
                        .[.[["nGene"]] >= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            metrics <- metrics %>%
                .[.[["nGene"]] >= minGenes, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `minGenes` cutoff")
        }
        summaryCells[["minGenes"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("minGenes", ">=", min(minGenes)),
            sep = " | "
        )

        # maxGenes -------------------------------------------------------------
        if (!is.null(names(maxGenes))) {
            assert_are_set_equal(names(maxGenes), sampleNames)
            message(paste(
                "maxGenes: per sample mode",
                printString(maxGenes),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(maxGenes),
                cutoff = maxGenes,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleName"]] == sample, , drop = FALSE] %>%
                        .[.[["nGene"]] <= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            metrics <- metrics %>%
                .[.[["nGene"]] <= maxGenes, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `maxGenes` cutoff")
        }
        summaryCells[["maxGenes"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("maxGenes", "<=", max(maxGenes)),
            sep = " | "
        )

        # minNovelty -----------------------------------------------------------
        if (!is.null(names(minNovelty))) {
            assert_are_set_equal(names(minNovelty), sampleNames)
            message(paste(
                "minNovelty: per sample mode",
                printString(minNovelty),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(minNovelty),
                cutoff = minNovelty,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleName"]] == sample, , drop = FALSE] %>%
                        .[.[["log10GenesPerUMI"]] >= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            metrics <- metrics %>%
                .[.[["log10GenesPerUMI"]] >= minNovelty, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `minNovelty` cutoff")
        }
        summaryCells[["minNovelty"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("minNovelty", "<=", min(minNovelty)),
            sep = " | "
        )

        # maxMitoRatio ---------------------------------------------------------
        if (!is.null(names(maxMitoRatio))) {
            assert_are_set_equal(names(maxMitoRatio), sampleNames)
            message(paste(
                "maxMitoRatio: per sample mode",
                printString(maxMitoRatio),
                sep = "\n"
            ))
            list <- mapply(
                sample = names(maxMitoRatio),
                cutoff = maxMitoRatio,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleName"]] == sample, , drop = FALSE] %>%
                        .[.[["mitoRatio"]] >= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            metrics <- metrics %>%
                .[.[["mitoRatio"]] <= maxMitoRatio, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `maxMitoRatio` cutoff")
        }
        summaryCells[["maxMitoRatio"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("maxMitoRatio", "<=", max(maxMitoRatio)),
            sep = " | "
        )

        cells <- sort(rownames(metrics))
        assert_is_subset(cells, colnames(object))

        # Filter low quality genes =============================================
        summaryGenes <- character()
        summaryGenes[["prefilter"]] <- paste(
            paste(.paddedCount(nrow(object)), "genes"),
            "prefilter",
            sep = " | "
        )
        if (minCellsPerGene > 0L) {
            counts <- assay(object)
            numCells <- rowSums(counts > 0L)
            genes <- names(numCells[which(numCells >= minCellsPerGene)])
        } else {
            genes <- sort(rownames(object))
        }
        if (!length(genes)) {
            stop("No genes passed `minCellsPerGene` cutoff")
        }
        summaryGenes[["minCellsPerGene"]] <- paste(
            paste(.paddedCount(length(genes)), "genes"),
            paste("minCellsPerGene", "<=", as.character(minCellsPerGene)),
            sep = " | "
        )

        # Summary ==============================================================
        printParams <- c(
            paste(">=", min(minUMIs), "UMIs per cell"),
            paste("<=", max(maxUMIs), "UMIs per cell"),
            paste(">=", min(minGenes), "genes per cell"),
            paste("<=", max(maxGenes), "genes per cell"),
            paste(">=", min(minNovelty), "novelty score"),
            paste("<=", max(maxMitoRatio), "mitochondrial abundance"),
            paste(">=", min(minCellsPerGene), "cells per gene")
        )
        cat(c(
            "Parameters:",
            paste("  -", printParams),
            bcbioBase::separatorBar,
            "Cells:",
            as.character(summaryCells),
            bcbioBase::separatorBar,
            "Genes:",
            as.character(summaryGenes),
            bcbioBase::separatorBar,
            "Summary:",
            paste(
                "  -",
                length(cells), "of", ncol(object), "cells passed filtering",
                paste0("(", percent(length(cells) / ncol(object)), ")")
            ),
            paste(
                "  -",
                length(genes), "of", nrow(object), "genes passed filtering",
                paste0("(", percent(length(genes) / nrow(object)), ")")
            )
        ), sep = "\n")

        summary <- list(
            "cells" = summaryCells,
            "genes" = summaryGenes
        )

        # Metadata =============================================================
        metadata(object)[["cellularBarcodes"]] <- NULL
        metadata(object)[["filterCells"]] <- cells
        metadata(object)[["filterGenes"]] <- genes
        metadata(object)[["filterParams"]] <- params
        metadata(object)[["filterSummary"]] <- summary

        object[genes, cells]
    }
)
