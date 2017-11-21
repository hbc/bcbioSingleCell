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
#' @param showReport Show summary statistics report and plots.
#' @param headerLevel RMarkdown header level, if `showReport = TRUE`.
#' @param drop Drop low quality cells from object after filtering.
#'   Enabled by default and generally recommended. If set `FALSE`, then the this
#'   function simply saves `filterCells`, `filterGenes`, and `filterParams` into
#'   the [metadata()] slot.
#'
#' @seealso [Seurat::CreateSeuratObject()].
#'
#' @return [bcbioSingleCell] object, with filtering information slotted into
#'   [metadata()] as `filterCells` and `filterParams`.
#'
#' @examples
#' \dontrun{
#' data(bcb)
#' bcb <- filterCells(bcb)
#' }
NULL



# Constructors ====
#' @importFrom basejump mdHeader mdList
#' @importFrom Matrix rowSums
#' @importFrom scales percent
#' @importFrom tibble rownames_to_column
.filterCells <- function(
    object,
    minUMIs = 1000,
    minGenes = 500,
    maxGenes = Inf,
    maxMitoRatio = 0.1,
    minNovelty = 0.8,
    minCellsPerGene = 3,
    showReport = FALSE,
    headerLevel = 2,
    drop = TRUE) {
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
    # Add support `nCount` filtering in a future update

    # Filter low quality cells ====
    # Don't use `subset()` here. That function uses non-standard evaluation that
    # should only be used interactively in a script.

    colData <- colData(object)
    message(paste(
        nrow(colData), "cells before filtering"
    ))

    # minUMIs
    if (!is.null(minUMIs) & minUMIs > 0) {
        colData <- colData %>%
            .[.[["nUMI"]] >= minUMIs, , drop = FALSE]
        message(paste(
            nrow(colData), "cells after 'minUMIs' cutoff"
        ))
    } else {
        message("'minUMIs' cutoff not applied")
    }

    # minGenes
    if (!is.null(minGenes) & minGenes > 0) {
        colData <- colData %>%
            .[.[["nGene"]] >= minGenes, , drop = FALSE]
        message(paste(
            nrow(colData), "cells after 'minGenes' cutoff"
        ))
    } else {
        message("'minGenes' cutoff not applied")
    }

    # maxGenes
    if (!is.null(maxGenes) & maxGenes < Inf) {
        colData <- colData %>%
            .[.[["nGene"]] <= maxGenes, , drop = FALSE]
        message(paste(
            nrow(colData), "cells after 'maxGenes' cutoff"
        ))
    } else {
        message("'maxGenes' cutoff not applied")
    }

    # maxMitoRatio
    if (!is.null(maxMitoRatio) & maxMitoRatio < 1) {
        colData <- colData %>%
            .[.[["mitoRatio"]] <= maxMitoRatio, , drop = FALSE]
        message(paste(
            nrow(colData), "cells after 'maxMitoRatio' cutoff"
        ))
    } else {
        message("'maxMitoRatio' cutoff not applied")
    }

    # minNovelty
    if (!is.null(minNovelty) & minNovelty > 0) {
        colData <- colData %>%
            .[.[["log10GenesPerUMI"]] >= minNovelty, , drop = FALSE]
        message(paste(
            nrow(colData), "cells after 'minNovelty' cutoff"
        ))
    } else {
        message("'minNovelty' cutoff not applied")
    }

    # Check for remaining cells
    if (!nrow(colData)) {
        stop("No cells passed filtering", call. = FALSE)
    }

    cells <- rownames(colData)
    message(paste(
        length(cells),
        "/",
        ncol(object),
        "cells passed filtering",
        paste0("(", percent(length(cells) / ncol(object)), ")")
    ))

    # Filter low expression genes ====
    if (minCellsPerGene > 0) {
        counts <- assay(object)
        numCells <- Matrix::rowSums(counts > 0)
        genes <- names(numCells[which(numCells >= minCellsPerGene)])
        message(paste(
            length(genes),
            "/",
            nrow(object),
            "genes passed filtering",
            paste0("(", percent(length(genes) / nrow(object)), ")")
        ))
    } else {
        genes <- NULL
    }

    # Metadata ====
    metadata(object)[["filterCells"]] <- cells
    metadata(object)[["filterGenes"]] <- genes
    metadata(object)[["filterParams"]] <- params
    cell2sample <- cell2sample(
        cells,
        samples = sampleMetadata(object)[["sampleID"]]
    )
    metadata(object)[["cell2sample"]] <- cell2sample

    # Drop cells and genes (destructive) ====
    if (isTRUE(drop)) {
        message(paste(
            "Dropping low quality cells and genes from the object"
        ))
        object <- .applyFilterCutoffs(object)
    } else {
        message(paste(
            "Non-destructive mode",
            "  cutoffs: metadata(object)$filterParams",
            "    cells: metadata(object)$filterCells",
            "    genes: metadata(object)$filterGenes",
            sep = "\n"
        ))
    }

    # Show summary statistics report and plots, if desired
    if (isTRUE(showReport)) {
        mdHeader("Filter parameters", level = headerLevel, asis = TRUE)
        mdList(c(
            paste0("`>= ", minUMIs, "` UMI counts per cell"),
            paste0("`>= ", minGenes, "` genes per cell"),
            paste0("`<= ", maxGenes, "` genes per cell"),
            paste0("`<= ", maxMitoRatio, "` relative mitochondrial abundance"),
            paste0("`>= ", minNovelty, "` novelty score")),
            asis = TRUE)

        mdHeader(
            "Filtered metrics plots",
            level = headerLevel,
            tabset = TRUE,
            asis = TRUE)

        # Reads per cell currently only supported for bcbio runs
        if (metadata(object)[["pipeline"]] == "bcbio") {
            mdHeader(
                "Reads per cell",
                level = headerLevel + 1,
                asis = TRUE)
            plotReadsPerCell(object) %>%
                show()
        }

        mdHeader(
            "Cell counts",
            level = headerLevel + 1,
            asis = TRUE)
        plotCellCounts(object) %>%
            show()

        mdHeader(
            "UMI counts per cell",
            level = headerLevel + 1,
            asis = TRUE)
        plotUMIsPerCell(object) %>%
            show()

        mdHeader(
            "Genes detected",
            level = headerLevel + 1,
            asis = TRUE)
        plotGenesPerCell(object) %>%
            show()

        mdHeader(
            "UMIs vs. genes",
            level = headerLevel + 1,
            asis = TRUE)
        plotUMIsVsGenes(object) %>%
            show()

        mdHeader(
            "Mitochondrial counts ratio",
            level = headerLevel + 1,
            asis = TRUE)
        plotMitoRatio(object) %>%
            show()

        mdHeader(
            "Novelty",
            level = headerLevel + 1,
            asis = TRUE)
        plotNovelty(object) %>%
            show()
    }

    object
}



# Methods ====
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
