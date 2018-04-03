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
#' @param minUMI Minimum number of UMI disambiguated counts per cell.
#' @param maxUMI Maximum number of UMI disambiguated counts per cell.
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
#' print(bcb_small)
#' filterCells(bcb_small, minGenes = 0L)
NULL



# Constructors =================================================================
.filterCells <- function(
    object,
    minUMI = 1000L,
    maxUMI = Inf,
    minGenes = 500L,
    maxGenes = Inf,
    maxMitoRatio = 0.1,
    minNovelty = 0.75,
    minCellsPerGene = 3L
) {
    validObject(object)
    metrics <- metrics(object)
    samples <- levels(metrics[["sampleID"]])

    # Legacy arguments =========================================================
    call <- match.call(expand.dots = TRUE)
    # maxUMIs
    if ("maxUMIs" %in% names(call)) {
        warn("Use `maxUMI` instead of `maxUMIs`")
        maxUMI <- call[["maxUMIs"]]
    }
    # minUMIs
    if ("minUMIs" %in% names(call)) {
        warn("Use `minUMI` instead of `minUMIs`")
        minUMI <- call[["minUMIs"]]
    }

    # Parameter integrity checks ===============================================
    params <- list(
        minUMI = minUMI,
        maxUMI = maxUMI,
        minGenes = minGenes,
        maxGenes = maxGenes,
        maxMitoRatio = maxMitoRatio,
        minNovelty = minNovelty,
        minCellsPerGene = minCellsPerGene
    )
    invisible(lapply(params, assert_is_a_number))
    assert_all_are_non_negative(as.numeric(params))

    # Filter low quality cells =================================================
    summary <- list()
    summary[["prefilter"]] <- paste(
        nrow(object), "genes",
        "/",
        ncol(object), "cells"
    )

    # minUMI -------------------------------------------------------------------
    if (!is.null(names(minUMI))) {
        # Per sample mode
        assert_is_subset(names(minUMI), samples)
        list <- mapply(
            sample = names(minUMI),
            cutoff = minUMI,
            FUN = function(sample, cutoff) {
                metrics %>%
                    .[.[["sampleID"]] == sample, , drop = FALSE] %>%
                    .[.[["nUMI"]] >= cutoff, , drop = FALSE]
            },
            SIMPLIFY = FALSE
        )
        metrics <- do.call(rbind, list)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nUMI"]] >= minUMI, , drop = FALSE]
    }
    assert_has_rows(metrics)
    summary[["minUMI"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("minUMI", ">=", toString(minUMI))
    )

    # maxUMI -------------------------------------------------------------------
    if (!is.null(names(maxUMI))) {
        # Per sample mode
        assert_is_subset(names(maxUMI), samples)
        list <- mapply(
            sample = names(maxUMI),
            cutoff = maxUMI,
            FUN = function(sample, cutoff) {
                metrics %>%
                    .[.[["sampleID"]] == sample, , drop = FALSE] %>%
                    .[.[["nUMI"]] <= cutoff, , drop = FALSE]
            },
            SIMPLIFY = FALSE
        )
        metrics <- do.call(rbind, list)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nUMI"]] <= maxUMI, , drop = FALSE]
    }
    assert_has_rows(metrics)
    summary[["maxUMI"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("maxUMI", "<=", toString(maxUMI))
    )

    # minGenes -----------------------------------------------------------------
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
            SIMPLIFY = FALSE
        )
        metrics <- do.call(rbind, list)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nGene"]] >= minGenes, , drop = FALSE]
    }
    assert_has_rows(metrics)
    summary[["minGenes"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("minGenes", ">=", toString(minGenes))
    )

    # maxGenes -----------------------------------------------------------------
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
            SIMPLIFY = FALSE
        )
        metrics <- do.call(rbind, list)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nGene"]] <= maxGenes, , drop = FALSE]
    }
    assert_has_rows(metrics)
    summary[["maxGenes"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("maxGenes", "<=", toString(maxGenes))
    )

    # maxMitoRatio -------------------------------------------------------------
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
            SIMPLIFY = FALSE
        )
        metrics <- do.call(rbind, list)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["mitoRatio"]] <= maxMitoRatio, , drop = FALSE]
    }
    assert_has_rows(metrics)
    summary[["maxMitoRatio"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("maxMitoRatio", "<=", toString(maxMitoRatio))
    )

    # minNovelty ---------------------------------------------------------------
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
            SIMPLIFY = FALSE
        )
        metrics <- do.call(rbind, list)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["log10GenesPerUMI"]] >= minNovelty, , drop = FALSE]
    }
    assert_has_rows(metrics)
    summary[["minNovelty"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("minNovelty", "<=", toString(minNovelty))
    )

    cells <- sort(rownames(metrics))
    assert_is_subset(cells, colnames(object))

    # Filter low quality genes =================================================
    if (minCellsPerGene > 0L) {
        counts <- assay(object)
        numCells <- Matrix::rowSums(counts > 0L)
        genes <- names(numCells[which(numCells >= minCellsPerGene)])
    } else {
        genes <- sort(rownames(object))
    }
    assert_is_non_empty(genes)
    summary[["minCellsPerGene"]] <- paste(
        paste(.paddedCount(length(genes)), "genes"),
        "|",
        paste("minCellsPerGene", "<=", as.character(minCellsPerGene))
    )

    # Summary ==================================================================
    printParams <- c(
        paste(">=", toString(minUMI), "UMI per cell"),
        paste("<=", toString(maxUMI), "UMI per cell"),
        paste(">=", toString(minGenes), "genes per cell"),
        paste("<=", toString(maxGenes), "genes per cell"),
        paste("<=", toString(maxMitoRatio), "mitochondrial abundance"),
        paste(">=", toString(minNovelty), "novelty score"),
        paste(">=", toString(minCellsPerGene), "cells per gene")
    )
    cat(c(
        "Filtering parameters:",
        paste("  -", printParams),
        bcbioBase::separatorBar,
        as.character(summary),
        bcbioBase::separatorBar,
        paste(
            length(cells), "/", ncol(object), "cells passed filtering",
            paste0("(", percent(length(cells) / ncol(object)), ")")
        ),
        paste(
            length(genes), "/", nrow(object), "genes passed filtering",
            paste0("(", percent(length(genes) / nrow(object)), ")")
        )
    ), sep = "\n")

    # Metadata =================================================================
    metadata(object)[["cellularBarcodes"]] <- NULL
    metadata(object)[["filterCells"]] <- cells
    metadata(object)[["filterGenes"]] <- genes
    metadata(object)[["filterParams"]] <- params
    metadata(object)[["filterSummary"]] <- summary
    metadata <- Filter(Negate(is.null), metadata)

    .applyFilterCutoffs(object)
}



.paddedCount <- function(x, width = 8L) {
    str_pad(x, width = width, pad = " ")
}



# Methods ======================================================================
#' @rdname filterCells
#' @export
setMethod(
    "filterCells",
    signature("bcbioSingleCell"),
    .filterCells
)
