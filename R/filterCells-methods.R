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
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param nCells `scalar integer`. Expected number of cells per sample.
#' @param minUMIs `scalar integer`. Minimum number of UMI disambiguated counts
#'   per cell.
#' @param maxUMIs `scalar integer`. Maximum number of UMI disambiguated counts
#'   per cell.
#' @param minGenes `scalar integer`. Minimum number of genes detected.
#' @param maxGenes `scalar integer`. Maximum number of genes detected.
#' @param minNovelty `scalar integer` (`0`-`1`). Minimum novelty score (log10
#'   genes per UMI).
#' @param maxMitoRatio `scalar integer` (`0`-`1`). Maximum relative
#'   mitochondrial abundance.
#' @param minCellsPerGene `scalar integer`. Include genes with non-zero
#'   expression in at least this many cells.
#'
#' @seealso [Seurat::CreateSeuratObject()].
#'
#' @return `bcbioSingleCell`, with filtering information slotted into
#'   [metadata()] as `filterCells` and `filterParams`.
#'
#' @examples
#' # SingleCellExperiment ====
#' object <- cellranger_small
#' show(object)
#'
#' x <- filterCells(object)
#' show(x)
#' metadata(x)$filterParams
#'
#' # Per sample cutoffs
#' sampleNames(object)
#' x <- filterCells(
#'     object = object,
#'     minUMIs = c(pbmc = 100)
#' )
#' show(x)
#' metadata(x)$filterParams
NULL



# Constructors =================================================================
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



# Methods ======================================================================
#' @rdname filterCells
#' @export
setMethod(
    "filterCells",
    signature("SingleCellExperiment"),
    function(
        object,
        nCells = Inf,
        minUMIs = 0L,
        maxUMIs = Inf,
        minGenes = 0L,
        maxGenes = Inf,
        minNovelty = 0L,
        maxMitoRatio = 1L,
        minCellsPerGene = 1L,
        zinbwave = FALSE
    ) {
        validObject(object)
        assert_is_a_bool(zinbwave)

        sampleNames <- sampleNames(object)
        metrics <- metrics(object)

        # Parameter integrity checks ===========================================
        # Expected nCells per sample
        assert_is_a_number(nCells)
        assert_all_are_positive(nCells)

        # minUMIs
        assert_is_any_of(minUMIs, c("numeric", "character"))
        if (is.character(minUMIs)) {
            assert_is_a_string(minUMIs)
            assert_is_subset(minUMIs, c("inflection", "knee"))
        }

        # maxUMIs
        assert_is_numeric(maxUMIs)
        assert_all_are_positive(maxUMIs)

        # minGenes
        assert_is_numeric(minGenes)
        assert_all_are_non_negative(minGenes)

        # maxGenes
        assert_is_numeric(maxGenes)
        assert_all_are_non_negative(maxGenes)

        # minNovelty
        assert_is_numeric(minNovelty)
        assert_all_are_in_range(minNovelty, lower = 0L, upper = 1L)

        # maxMitoRatio
        assert_is_numeric(maxMitoRatio)
        assert_all_are_in_range(maxMitoRatio, lower = 0L, upper = 1L)

        # minCellsPerGene
        # Don't allow genes with all zero counts, so require at least 1 here
        assert_is_numeric(minCellsPerGene)
        assert_all_are_positive(minCellsPerGene)

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

        # Expected nCells per sample (filtered by top nUMI) --------------------
        if (nCells < Inf) {
            metrics <- metrics %>%
                rownames_to_column() %>%
                group_by(!!sym("sampleID")) %>%
                arrange(desc(!!sym("nUMI")), .by_group = TRUE) %>%
                slice(seq_len(nCells)) %>%
                as.data.frame() %>%
                column_to_rownames()
        }

        if (!nrow(metrics)) {
            stop("No cells passed `nCells` cutoff")
        }
        summaryCells[["nCells"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("nCells", "==", nCells),
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
            nonzeroGenes <- rowSums(counts(object) > 0L)
            keep <- nonzeroGenes[which(nonzeroGenes >= minCellsPerGene)]
            genes <- names(keep)
        } else {
            genes <- rownames(object)
        }
        genes <- sort(genes)
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
            paste("==", nCells, "cells per sample"),
            paste(">=", min(minCellsPerGene), "cells per gene")
        )
        message(paste(c(
            "Parameters:",
            paste("  -", printParams),
            bcbioBase::separatorBar,
            "Cells:",
            as.character(summaryCells),
            printString(table(metrics[["sampleName"]])),
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
        ), collapse = "\n"))

        summary <- list(cells = summaryCells, genes = summaryGenes)

        # Metadata =============================================================
        metadata(object)[["cellularBarcodes"]] <- NULL
        metadata(object)[["filterCells"]] <- cells
        metadata(object)[["filterGenes"]] <- genes
        metadata(object)[["filterParams"]] <- params
        metadata(object)[["filterSummary"]] <- summary

        # Subset ===============================================================
        object <- object[genes, cells]

        # zinbwave weights =====================================================
        if (isTRUE(zinbwave)) {
            object <- .slotZinbwaveIntoAssays(object)
        }

        object
    }
)
