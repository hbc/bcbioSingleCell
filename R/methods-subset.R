#' Bracket-Based Subsetting
#'
#' Extract genes by row and cells by column from a `bcbioSingleCell` object.
#'
#' @name subset
#' @author Michael Steinbaugh
#'
#' @importFrom parallel mclapply
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
#' cells <- head(colnames(bcb), 100L)
#' head(cells)
#' genes <- head(rownames(bcb), 100L)
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
.subset <- function(x, i, j, ..., drop = FALSE) {
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

    assays <- assays(sce)
    rowRanges <- rowRanges(sce)
    colData <- colData(sce)
    isSpike <- isSpike(sce)

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

    # Unfiltered cellular barcodes
    cb <- metadata[["cellularBarcodes"]]
    # Bind barcodes into a single `data.frame`, which we can subset
    if (!is.null(cb)) {
        assert_is_list(cb)
        df <- cb %>%
            .bindCellularBarcodes() %>%
            .[cells, , drop = FALSE]
        cb <- mclapply(seq_along(sampleIDs), function(a) {
            df %>%
                ungroup() %>%
                filter(.data[["sampleID"]] == sampleIDs[[a]]) %>%
                mutate(sampleID = NULL)
        })
        names(cb) <- sampleIDs
        metadata[["cellularBarcodes"]] <- cb
    }

    # Return ===================================================================
    sce <- SingleCellExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata
    )
    # Define spikeNames for spike-in sequences
    if (is.character(isSpike)) {
        for (i in seq_along(isSpike)) {
            isSpike(sce, isSpike[[i]]) <- isSpike[[i]]
        }
    }
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
    .subset
)
