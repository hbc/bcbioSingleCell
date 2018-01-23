#' Filter Cells
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @details The filtering cutoff values now support a named numeric vector. By
#' default we recommend applying the same filtering cutoff to all samples.
#' When matching the samples, be sure to use `sampleID` column (i.e. the
#' rownames of [sampleMetadata()]).
#'
#' @rdname filterCells
#' @name filterCells
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param minUMIs Minimum number of UMI disambiguated counts per cell.
#' @param maxUMIs Maximum number of UMI disambiguated counts per cell.
#' @param minGenes Minimum number of genes detected.
#' @param maxGenes Maximum number of genes detected.
#' @param maxMitoRatio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param minNovelty Minimum novelty score.
#' @param minCellsPerGene Include genes with non-zero expression in at least
#'   this many cells.
#' @param quiet If `TRUE`, don't show the filtering parameter summary.
#'
#' @seealso [Seurat::CreateSeuratObject()].
#'
#' @return [bcbioSingleCell] object, with filtering information slotted into
#'   [metadata()] as `filterCells` and `filterParams`.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' print(bcb)
#' filterCells(bcb)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase mdHeader mdList
#' @importFrom Matrix rowSums
#' @importFrom scales percent
#' @importFrom tibble rownames_to_column
.filterCells <- function(
    object,
    minUMIs = 1000,
    maxUMIs = Inf,
    minGenes = 500,
    maxGenes = Inf,
    maxMitoRatio = 0.1,
    minNovelty = 0.75,
    minCellsPerGene = 3,
    quiet = FALSE) {
    metrics <- metrics(object)

    # Parameter integrity checks ===============================================
    params <- list(
        minUMIs = minUMIs,
        maxUMIs = maxUMIs,
        minGenes = minGenes,
        maxGenes = maxGenes,
        maxMitoRatio = maxMitoRatio,
        minNovelty = minNovelty,
        minCellsPerGene = minCellsPerGene)
    # Ensure all params are numeric
    if (!all(vapply(
        X = params,
        FUN = is.numeric,
        FUN.VALUE = logical(1L)
    ))) {
        abort("Filter parameters must be numeric")
    }
    # Ensure all params are not negative
    if (!all(vapply(
        X = params,
        FUN = function(x) { all(x >= 0L) },
        FUN.VALUE = logical(1L)
    ))) {
        abort("Filter parameters must be non-negative")
    }

    # Filter low quality cells =================================================
    summary <- list()
    summary[["prefilter"]] <- paste(
        nrow(object), "genes",
        "/",
        ncol(object), "cells"
    )


    # minUMIs ====
    if (!is.null(names(minUMIs))) {
        # Per sample mode
        if (!all(names(minUMIs) %in% metrics[["sampleID"]])) {
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
    summary[["minUMIs"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("minUMIs", ">=", toString(minUMIs))
    )

    # maxUMIs ====
    if (!is.null(names(maxUMIs))) {
        # Per sample mode
        if (!all(names(maxUMIs) %in% metrics[["sampleID"]])) {
            abort("`maxUMIs` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(maxUMIs)), function(a) {
            sampleID <- names(maxUMIs)[[a]]
            cutoff <- maxUMIs[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["nUMI"]] <= cutoff, , drop = FALSE]
        })
        metrics <- do.call(rbind, metrics)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nUMI"]] <= maxUMIs, , drop = FALSE]
    }
    summary[["maxUMIs"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("maxUMIs", "<=", toString(maxUMIs))
    )

    # minGenes ====
    if (!is.null(names(minGenes))) {
        # Per sample mode
        if (!all(names(minGenes) %in% metrics[["sampleID"]])) {
            abort("`minGenes` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(minGenes)), function(a) {
            sampleID <- names(minGenes)[[a]]
            cutoff <- minGenes[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["nGene"]] >= cutoff, , drop = FALSE]
        })
        metrics <- do.call(rbind, metrics)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nGene"]] >= minGenes, , drop = FALSE]
    }
    summary[["minGenes"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("minGenes", ">=", toString(minGenes))
    )

    # maxGenes ====
    if (!is.null(names(maxGenes))) {
        # Per sample mode
        if (!all(names(maxGenes) %in% metrics[["sampleID"]])) {
            abort("`maxGenes` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(maxGenes)), function(a) {
            sampleID <- names(maxGenes)[[a]]
            cutoff <- maxGenes[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["nGene"]] <= cutoff, , drop = FALSE]
        })
        metrics <- do.call(rbind, metrics)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["nGene"]] <= maxGenes, , drop = FALSE]
    }
    summary[["maxGenes"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("maxGenes", "<=", toString(maxGenes))
    )

    # maxMitoRatio ====
    if (!is.null(names(maxMitoRatio))) {
        # Per sample mode
        if (!all(names(maxMitoRatio) %in% metrics[["sampleID"]])) {
            abort("`maxMitoRatio` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(maxMitoRatio)), function(a) {
            sampleID <- names(maxMitoRatio)[[a]]
            cutoff <- maxMitoRatio[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["mitoRatio"]] >= cutoff, , drop = FALSE]
        })
        metrics <- do.call(rbind, metrics)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["mitoRatio"]] <= maxMitoRatio, , drop = FALSE]
    }
    summary[["maxMitoRatio"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("maxMitoRatio", "<=", toString(maxMitoRatio))
    )

    # minNovelty ====
    if (!is.null(names(minNovelty))) {
        # Per sample mode
        if (!all(names(minNovelty) %in% metrics[["sampleID"]])) {
            abort("`minNovelty` names don't match sample IDs")
        }
        list <- lapply(seq_along(names(minNovelty)), function(a) {
            sampleID <- names(minNovelty)[[a]]
            cutoff <- minNovelty[[a]]
            metrics %>%
                .[.[["sampleID"]] == sampleID, , drop = FALSE] %>%
                .[.[["log10GenesPerUMI"]] >= cutoff, , drop = FALSE]
        })
        metrics <- do.call(rbind, metrics)
    } else {
        # Fixed cutoff value
        metrics <- metrics %>%
            .[.[["log10GenesPerUMI"]] >= minNovelty, , drop = FALSE]
    }
    summary[["minNovelty"]] <- paste(
        paste(.paddedCount(nrow(metrics)), "cells"),
        "|",
        paste("minNovelty", "<=", toString(minNovelty))
    )

    cells <- sort(rownames(metrics))
    if (!length(cells)) {
        warn("No cells passed filtering cutoffs")
        return(NULL)
    }
    # Check to make sure the cells are valid
    if (!all(cells %in% colnames(object))) {
        # The tidyverse chain of tools has a tendency to drop rownames. Be
        # sure to use base R methods above for our filtering cutoffs.
        abort("Cell vector unexpectedly doens't match IDs in object colnames")
    }

    # Filter low quality genes =================================================
    if (minCellsPerGene > 0) {
        counts <- assay(object)
        numCells <- Matrix::rowSums(counts > 0)
        genes <- names(numCells[which(numCells >= minCellsPerGene)])
    } else {
        genes <- sort(rownames(object))
    }
    summary[["minCellsPerGene"]] <- paste(
        paste(.paddedCount(length(genes)), "genes"),
        "|",
        paste("minCellsPerGene", "<=", as.character(minCellsPerGene))
    )
    if (!length(genes)) {
        warn("No genes passed filtering")
        return(NULL)
    }

    # Summary ==================================================================
    if (!isTRUE(quiet)) {
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
            sepBar,
            as.character(summary),
            sepBar,
            paste(
                length(cells), "/", ncol(object), "cells passed filtering",
                paste0("(", percent(length(cells) / ncol(object)), ")")
            ),
            paste(
                length(genes), "/", nrow(object), "genes passed filtering",
                paste0("(", percent(length(genes) / nrow(object)), ")")
            )
        ), sep = "\n")
    }

    # Metadata =================================================================
    metadata(object)[["filterCells"]] <- cells
    metadata(object)[["filterGenes"]] <- genes
    metadata(object)[["filterParams"]] <- params
    metadata(object)[["filterSummary"]] <- summary

    .applyFilterCutoffs(object)
}



#' @importFrom stringr str_pad
.paddedCount <- function(x, width = 8) {
    str_pad(x, width = width, pad = " ")
}



# Methods ======================================================================
#' @rdname filterCells
#' @export
setMethod(
    "filterCells",
    signature("bcbioSingleCell"),
    .filterCells)



#' @rdname filterCells
#' @export
setMethod(
    "filterCells",
    signature("bcbioSCDataSet"),
    function(object) {
        abort(paste(
            "Convert 'bcbioSCDataSet' to 'bcbioSingleCell' class.\n",
            "Run this code: bcb <- as(bcb, \"bcbioSingleCell\")"
        ))
    })
