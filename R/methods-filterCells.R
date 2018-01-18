#' Filter Cells
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filterCells
#' @name filterCells
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param minUMIs Minimum number of UMI disambiguated counts per cell.
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
    minGenes = 500,
    maxGenes = Inf,
    maxMitoRatio = 0.1,
    minNovelty = 0.75,
    minCellsPerGene = 3,
    quiet = FALSE) {
    # Ensure that all filter parameters are numeric
    params <- c(
        minUMIs = minUMIs,
        minGenes = minGenes,
        maxGenes = maxGenes,
        maxMitoRatio = maxMitoRatio,
        minNovelty = minNovelty,
        minCellsPerGene = minCellsPerGene)
    if (!is.numeric(params)) {
        stop("Filter parameters must be numeric", call. = FALSE)
    }
    # TODO Add support `nCount` filtering in a future update

    # Filter low quality cells =================================================
    colData <- colData(object)
    if (!isTRUE(quiet)) {
        message(paste(
            paste(ncol(object), "cells before filtering"),
            paste(nrow(object), "genes before filtering"),
            sep = "\n"
        ))
    }

    # minUMIs
    if (!is.null(minUMIs) & minUMIs > 0) {
        colData <- colData %>%
            .[.[["nUMI"]] >= minUMIs, , drop = FALSE]
    }
    if (!isTRUE(quiet)) {
        message(paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            "|",
            paste("minUMIs", ">=", as.character(minUMIs))
        ))
    }

    # minGenes
    if (!is.null(minGenes) & minGenes > 0) {
        colData <- colData %>%
            .[.[["nGene"]] >= minGenes, , drop = FALSE]
    }
    if (!isTRUE(quiet)) {
        message(paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            "|",
            paste("minGenes", ">=", as.character(minGenes))
        ))
    }

    # maxGenes
    if (!is.null(maxGenes) & maxGenes < Inf) {
        colData <- colData %>%
            .[.[["nGene"]] <= maxGenes, , drop = FALSE]
    }
    if (!isTRUE(quiet)) {
        message(paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            "|",
            paste("maxGenes", "<=", as.character(maxGenes))
        ))
    }

    # maxMitoRatio
    if (!is.null(maxMitoRatio) & maxMitoRatio < 1) {
        colData <- colData %>%
            .[.[["mitoRatio"]] <= maxMitoRatio, , drop = FALSE]
    }
    if (!isTRUE(quiet)) {
        message(paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            "|",
            paste("maxMitoRatio", "<=", as.character(maxMitoRatio))
        ))
    }

    # minNovelty
    if (!is.null(minNovelty) & minNovelty > 0) {
        colData <- colData %>%
            .[.[["log10GenesPerUMI"]] >= minNovelty, , drop = FALSE]
    }
    if (!isTRUE(quiet)) {
        message(paste(
            paste(.paddedCount(nrow(colData)), "cells"),
            "|",
            paste("minNovelty", "<=", as.character(minNovelty))
        ))
    }

    cells <- rownames(colData)
    if (!length(cells)) {
        warning("No cells passed filtering", call. = FALSE)
        return(NULL)
    }

    # Filter low quality genes =================================================
    if (minCellsPerGene > 0) {
        counts <- assay(object)
        numCells <- Matrix::rowSums(counts > 0)
        genes <- names(numCells[which(numCells >= minCellsPerGene)])
    } else {
        genes <- rownames(object)
    }
    if (!isTRUE(quiet)) {
        message(paste(
            paste(.paddedCount(length(genes)), "genes"),
            "|",
            paste("minCellsPerGene", "<=", as.character(minCellsPerGene))
        ))
    }
    if (!length(genes)) {
        warning("No genes passed filtering", call. = FALSE)
        return(NULL)
    }

    # Metadata =================================================================
    metadata(object)[["filterCells"]] <- cells
    metadata(object)[["filterGenes"]] <- genes
    metadata(object)[["filterParams"]] <- params

    # Summary ==================================================================
    if (!isTRUE(quiet)) {
        message(paste(
            paste(
                length(cells), "/", ncol(object), "cells passed filtering",
                paste0("(", percent(length(cells) / ncol(object)), ")")
            ),
            paste(
                length(genes), "/", nrow(object), "genes passed filtering",
                paste0("(", percent(length(genes) / nrow(object)), ")")
            ),
            sep = "\n"
        ))
        c(paste(">=", minUMIs, "UMI counts per cell"),
          paste(">=", minGenes, "genes per cell"),
          paste("<=", maxGenes, "genes per cell"),
          paste("<=", maxMitoRatio, "mitochondrial abundance"),
          paste(">=", minNovelty, "novelty score"),
          paste(">=", minCellsPerGene, "cells per gene")) %>%
            paste("  -", .) %>%
            c("Filtering parameters:", .) %>%
            cat(sep = "\n")
    }

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
        stop(paste(
            "Convert 'bcbioSCDataSet' to 'bcbioSingleCell' class.\n",
            "Run this code: bcb <- as(bcb, \"bcbioSingleCell\")"
        ), call. = FALSE)
})
