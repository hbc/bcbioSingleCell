#' Bracket-Based Subsetting
#'
#' Extract genes by row and cells by column from a `bcbioSingleCell` object.
#'
#' @name subset
#' @author Michael Steinbaugh
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @inheritParams base::`[`
#'
#' @param ... Additional arguments.
#'
#' @seealso
#' - `help("[", "base")`.
#' - `selectSamples()` for subsetting based on sample metadata.
#'
#' @return `bcbioSingleCell`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#'
#' cells <- colnames(bcb)[1:100]
#' head(cells)
#' genes <- rownames(bcb)[1:100]
#' head(genes)
#'
#' # Subset by cell identifiers
#' bcb[, cells]
#'
#' # Subset by genes
#' bcb[genes, ]
#'
#' # Subset by both genes and cells
#' bcb[genes, cells]
NULL



# Constructors =================================================================
#' @importFrom dplyr mutate_all
#' @importFrom magrittr set_rownames
#' @importFrom S4Vectors metadata SimpleList
.subset <- function(x, i, j, ..., drop) {
    validObject(x)

    # Genes
    if (missing(i)) {
        i <- 1L:nrow(x)
    }
    # Cells
    if (missing(j)) {
        j <- 1L:ncol(x)
    }

    # Early return if dimensions are unmodified
    if (identical(dim(x), c(length(i), length(j)))) {
        return(x)
    }

    # Regenerate and subset SummarizedExperiment
    sce <- as(x, "SingleCellExperiment")
    sce <- sce[i, j, drop = drop]

    genes <- rownames(sce)
    cells <- colnames(sce)

    # Assays ===================================================================
    assays <- assays(sce)

    # Row data =================================================================
    rowRanges <- rowRanges(sce)
    assert_has_names(rowRanges)

    # Column data ==============================================================
    # Don't need to relevel factors here currently
    colData <- colData(sce)

    # Metadata =================================================================
    metadata <- metadata(sce)
    metadata[["subset"]] <- TRUE
    # Update version, if necessary
    if (!identical(metadata[["version"]], packageVersion)) {
        metadata[["originalVersion"]] <- metadata[["version"]]
        metadata[["version"]] <- packageVersion
    }

    # cell2sample
    cell2sample <- metadata[["cell2sample"]]
    assert_is_factor(cell2sample)
    # Note that we're subsetting `sampleData` by the factor levels in
    # `cell2sample`, so this must come first
    cell2sample <- droplevels(cell2sample[cells])
    metadata[["cell2sample"]] <- cell2sample

    # sampleData
    sampleData <- sampleData(x) %>%
        .[levels(cell2sample), , drop = FALSE] %>%
        mutate_all(as.factor) %>%
        mutate_all(droplevels) %>%
        set_rownames(.[["sampleID"]])
    metadata[["sampleData"]] <- sampleData
    sampleIDs <- as.character(sampleData[["sampleID"]])

    # aggregateReplicates
    aggregateReplicates <- metadata[["aggregateReplicates"]]
    if (!is.null(aggregateReplicates)) {
        intersect <- intersect(cells, aggregateReplicates)
        aggregateReplicates <- aggregateReplicates %>%
            .[. %in% intersect]
        metadata[["aggregateReplicates"]] <- aggregateReplicates
    }

    # filterCells
    filterCells <- metadata[["filterCells"]]
    if (!is.null(filterCells)) {
        filterCells <- intersect(filterCells, cells)
        metadata[["filterCells"]] <- filterCells
    }

    # filterGenes
    filterGenes <- metadata[["filterGenes"]]
    if (!is.null(filterGenes)) {
        filterGenes <- intersect(filterGenes, genes)
        metadata[["filterGenes"]] <- filterGenes
    }

    # bcbio ====================================================================
    # FIXME Need to move this into metadata and rework
    bcbio <- bcbio(x)
    if (is(bcbio, "SimpleList") & length(bcbio)) {
        # Cellular barcodes
        cb <- bcbio[["cellularBarcodes"]]
        # Bind barcodes into a single `data.frame`, which we can subset
        if (!is.null(cb)) {
            if (is.list(cb)) {
                cb <- .bindCellularBarcodes(cb)
            }
            cb <- cb[cells, , drop = FALSE]
            cbList <- lapply(seq_along(sampleIDs), function(a) {
                cb %>%
                    ungroup() %>%
                    filter(.data[["sampleID"]] == sampleIDs[[a]]) %>%
                    mutate(sampleID = NULL)
            })
            names(cbList) <- sampleIDs
        }
        # Return the SimpleList
        bcbio <- list(
            cellularBarcodes = cbList) %>%
            as("SimpleList")
    }

    # Return ===================================================================
    sce <- SingleCellExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata
    )
    new("bcbioSingleCell", sce)
}



# Methods ======================================================================
#' @rdname subset
#' @export
setMethod(
    "[",
    signature(
        x = "bcbioSingleCell",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    function(x, i, j, ..., drop = FALSE) {
        .subset(x, i, j, ..., drop)
    })
