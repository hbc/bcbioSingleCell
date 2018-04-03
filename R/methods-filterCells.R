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
#' print(bcb_small)
#' filterCells(bcb_small, minGenes = 0L)
NULL



# Constructors =================================================================
.filterCells <- function(
    object,
    minUMIs = 1000L,
    maxUMIs = Inf,
    minGenes = 500L,
    maxGenes = Inf,
    maxMitoRatio = 0.1,
    minNovelty = 0.75,
    minCellsPerGene = 3L
) {
    metrics <- metrics(object)
    sampleIDs <- levels(metrics[["sampleID"]])

    # Parameter integrity checks ===============================================
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

    # Filter low quality cells =================================================
    summary <- list()
    summary[["prefilter"]] <- paste(
        nrow(object), "genes",
        "/",
        ncol(object), "cells"
    )

    # minUMIs ------------------------------------------------------------------
    if (!is.null(names(minUMIs))) {
        # Per sample mode
        if (!all(names(minUMIs) %in% sampleIDs)) {
            abort("`minUMIs` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(minUMIs)), function(a) {
            sampleID <- names(minUMIs)[[a]]
            cutoff <- minUMIs[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["nUMI"]] >= cutoff, , drop = FALSE]
        })
        metrics <- do.call(rbind, list)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nUMI"]] >= minUMIs, , drop = FALSE]
    }
    assert_has_rows(metrics)
    summary[["minUMIs"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("minUMIs", ">=", toString(minUMIs))
    )

    # maxUMIs ------------------------------------------------------------------
    if (!is.null(names(maxUMIs))) {
        # Per sample mode
        if (!all(names(maxUMIs) %in% sampleIDs)) {
            abort("`maxUMIs` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(maxUMIs)), function(a) {
            sampleID <- names(maxUMIs)[[a]]
            cutoff <- maxUMIs[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["nUMI"]] <= cutoff, , drop = FALSE]
        })
        metrics <- do.call(rbind, list)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nUMI"]] <= maxUMIs, , drop = FALSE]
    }
    assert_has_rows(metrics)
    summary[["maxUMIs"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("maxUMIs", "<=", toString(maxUMIs))
    )

    # minGenes -----------------------------------------------------------------
    if (!is.null(names(minGenes))) {
        # Per sample mode
        if (!all(names(minGenes) %in% sampleIDs)) {
            abort("`minGenes` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(minGenes)), function(a) {
            sampleID <- names(minGenes)[[a]]
            cutoff <- minGenes[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["nGene"]] >= cutoff, , drop = FALSE]
        })
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
        if (!all(names(maxGenes) %in% sampleIDs)) {
            abort("`maxGenes` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(maxGenes)), function(a) {
            sampleID <- names(maxGenes)[[a]]
            cutoff <- maxGenes[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["nGene"]] <= cutoff, , drop = FALSE]
        })
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
        if (!all(names(maxMitoRatio) %in% sampleIDs)) {
            abort("`maxMitoRatio` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(maxMitoRatio)), function(a) {
            sampleID <- names(maxMitoRatio)[[a]]
            cutoff <- maxMitoRatio[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["mitoRatio"]] >= cutoff, , drop = FALSE]
        })
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
        if (!all(names(minNovelty) %in% sampleIDs)) {
            abort("`minNovelty` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(minNovelty)), function(a) {
            sampleID <- names(minNovelty)[[a]]
            cutoff <- minNovelty[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["log10GenesPerUMI"]] >= cutoff, , drop = FALSE]
        })
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
        paste(">=", toString(minUMIs), "UMI counts per cell"),
        paste("<=", toString(maxUMIs), "UMI counts per cell"),
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
