#' Filter Cells
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @details The filtering cutoff values now support a named numeric vector. By
#' default we recommend applying the same filtering cutoff to all samples.
#' When matching the samples, be sure to use `sampleID` column (i.e. the
#' rownames of [sampleData()]).
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
#' @param maxMitoRatio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param minNovelty Minimum novelty score.
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
#' show(bcb_small)
#' x <- filterCells(bcb_small, minGenes = 0L)
#' show(x)
#'
#' # SingleCellExperiment ====
#' show(cellranger_small)
#' x <- filterCells(cellranger_small, minNovelty = 0L)
#' show(x)
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
        minUMIs = 1000L,
        maxUMIs = Inf,
        minGenes = 100L,
        maxGenes = Inf,
        maxMitoRatio = 0.1,
        minNovelty = 0.7,
        minCellsPerGene = 3L
    ) {
        validObject(object)
        metrics <- metrics(object)
        samples <- levels(metrics[["sampleID"]])

        # Parameter integrity checks ===========================================
        params <- list(
            minUMIs = minUMIs,
            maxUMIs = maxUMIs,
            minGenes = minGenes,
            maxGenes = maxGenes,
            maxMitoRatio = maxMitoRatio,
            minNovelty = minNovelty,
            minCellsPerGene = minCellsPerGene
        )
        invisible(lapply(params, assert_is_a_number))
        assert_all_are_non_negative(as.numeric(params))

        # Filter low quality cells =============================================
        summaryCells <- character()
        summaryCells[["prefilter"]] <- paste(
            paste(.paddedCount(ncol(object)), "cells"),
            "prefilter",
            sep = " | "
        )

        # minUMIs --------------------------------------------------------------
        if (!is.null(names(minUMIs))) {
            # Per sample mode
            assert_is_subset(names(minUMIs), samples)
            list <- mapply(
                sample = names(minUMIs),
                cutoff = minUMIs,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleID"]] == sample, , drop = FALSE] %>%
                        .[.[["nUMI"]] >= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            # Fixed cutoff value
            metrics <- metrics %>%
                .[.[["nUMI"]] >= minUMIs, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `minUMIs` cutoff")
        }
        summaryCells[["minUMIs"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("minUMIs", ">=", toString(minUMIs)),
            sep = " | "
        )

        # maxUMIs --------------------------------------------------------------
        if (!is.null(names(maxUMIs))) {
            # Per sample mode
            assert_is_subset(names(maxUMIs), samples)
            list <- mapply(
                sample = names(maxUMIs),
                cutoff = maxUMIs,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleID"]] == sample, , drop = FALSE] %>%
                        .[.[["nUMI"]] <= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            # Fixed cutoff value
            metrics <- metrics %>%
                .[.[["nUMI"]] <= maxUMIs, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `maxUMIs` cutoff")
        }
        summaryCells[["maxUMIs"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("maxUMIs", "<=", toString(maxUMIs)),
            sep = " | "
        )

        # minGenes -------------------------------------------------------------
        if (!is.null(names(minGenes))) {
            # Per sample mode
            assert_is_subset(names(minGenes), samples)
            list <- mapply(
                sample = names(minGenes),
                cutoff = minGenes,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleID"]] == sample, , drop = FALSE] %>%
                        .[.[["nGene"]] >= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            # Fixed cutoff value
            metrics <- metrics %>%
                .[.[["nGene"]] >= minGenes, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `minGenes` cutoff")
        }
        summaryCells[["minGenes"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("minGenes", ">=", toString(minGenes)),
            sep = " | "
        )

        # maxGenes -------------------------------------------------------------
        if (!is.null(names(maxGenes))) {
            # Per sample mode
            assert_is_subset(names(maxGenes), samples)
            list <- mapply(
                sample = names(maxGenes),
                cutoff = maxGenes,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleID"]] == sample, , drop = FALSE] %>%
                        .[.[["nGene"]] <= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            # Fixed cutoff value
            metrics <- metrics %>%
                .[.[["nGene"]] <= maxGenes, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `maxGenes` cutoff")
        }
        summaryCells[["maxGenes"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("maxGenes", "<=", toString(maxGenes)),
            sep = " | "
        )

        # maxMitoRatio ---------------------------------------------------------
        if (!is.null(names(maxMitoRatio))) {
            # Per sample mode
            assert_is_subset(names(maxMitoRatio), samples)
            list <- mapply(
                sample = names(maxMitoRatio),
                cutoff = maxMitoRatio,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleID"]] == sample, , drop = FALSE] %>%
                        .[.[["mitoRatio"]] >= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            # Fixed cutoff value
            metrics <- metrics %>%
                .[.[["mitoRatio"]] <= maxMitoRatio, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `maxMitoRatio` cutoff")
        }
        summaryCells[["maxMitoRatio"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("maxMitoRatio", "<=", toString(maxMitoRatio)),
            sep = " | "
        )

        # minNovelty -----------------------------------------------------------
        if (!is.null(names(minNovelty))) {
            # Per sample mode
            assert_is_subset(names(minNovelty), samples)
            list <- mapply(
                sample = names(minNovelty),
                cutoff = minNovelty,
                FUN = function(sample, cutoff) {
                    metrics %>%
                        .[.[["sampleID"]] == sample, , drop = FALSE] %>%
                        .[.[["log10GenesPerUMI"]] >= cutoff, , drop = FALSE]
                },
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )
            metrics <- do.call(rbind, list)
        } else {
            # Fixed cutoff value
            metrics <- metrics %>%
                .[.[["log10GenesPerUMI"]] >= minNovelty, , drop = FALSE]
        }
        if (!nrow(metrics)) {
            stop("No cells passed `minNovelty` cutoff")
        }
        summaryCells[["minNovelty"]] <- paste(
            paste(.paddedCount(nrow(metrics)), "cells"),
            paste("minNovelty", "<=", toString(minNovelty)),
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
            paste(">=", toString(minUMIs), "UMIs per cell"),
            paste("<=", toString(maxUMIs), "UMIs per cell"),
            paste(">=", toString(minGenes), "genes per cell"),
            paste("<=", toString(maxGenes), "genes per cell"),
            paste("<=", toString(maxMitoRatio), "mitochondrial abundance"),
            paste(">=", toString(minNovelty), "novelty score"),
            paste(">=", toString(minCellsPerGene), "cells per gene")
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

        .applyFilterCutoffs(object)
    }
)
