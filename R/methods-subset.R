#' Bracket-Based Subsetting
#'
#' Extract genes by row and cells by column from a [bcbioSingleCell] object.
#'
#' @rdname subset
#' @name subset
#'
#' @author Michael Steinbaugh
#'
#' @inheritParams base::`[`
#'
#' @param ... Additional arguments.
#'
#' @seealso
#' - `help("[", "base")`.
#' - `selectSamples()` for subsetting based on sample metadata.
#'
#' @return [bcbioSingleCell].
#'
#' @examples
#' # Subset by cell identifiers
#' \dontrun{
#' bcb[, cells]
#' }
#'
#' # Subset by genes
#' \dontrun{
#' bcb[genes, ]
#' }
#'
#' # Subset by both genes and cells
#' \dontrun{
#' bcb[genes, cells]
#' }
NULL



# Constructors ====
#' @importFrom dplyr mutate_if
#' @importFrom magrittr set_rownames
#' @importFrom S4Vectors metadata SimpleList
.subset <- function(x, i, j, ..., drop = FALSE) {
    if (missing(i)) {
        i <- 1:nrow(x)
    }
    if (missing(j)) {
        j <- 1:ncol(x)
    }

    # Prepare and subset SummarizedExperiment
    se <- as(x, "SummarizedExperiment")
    se <- se[i, j, drop = drop]

    genes <- rownames(se)
    cells <- colnames(se)

    assays <- assays(se)
    rowData <- rowData(se)
    if (!is.null(rowData)) {
        rownames(rowData) <- slot(se, "NAMES")
    }
    colData <- colData(se)

    # Metadata =================================================================
    metadata <- metadata(se)

    # cell2sample mappings
    cell2sample <- cell2sample(
        colnames(se),
        samples = rownames(sampleMetadata(x))
    )
    metadata[["cell2sample"]] <- cell2sample

    # sampleMetadata
    sampleMetadata <- sampleMetadata(x) %>%
        .[levels(cell2sample), ] %>%
        mutate_if(is.factor, droplevels) %>%
        set_rownames(.[["sampleID"]])
    metadata[["sampleMetadata"]] <- sampleMetadata
    sampleIDs <- as.character(sampleMetadata[["sampleID"]])

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

    metadata[["subset"]] <- TRUE

    # bcbio slot ===============================================================
    bcbio <- bcbio(x)
    if (!is.null(bcbio)) {
        # Cellular barcodes
        cb <- bcbio[["cellularBarcodes"]]
        # Bind barcodes into a single `data.frame`, which we can subset
        if (!is.null(cb)) {
            if (is.list(cb)) {
                cb <- .bindCellularBarcodes(cb)
            }
            cb <- cb[cells, ]
            cbList <- lapply(seq_along(sampleIDs), function(a) {
                cb %>%
                    ungroup() %>%
                    filter(sampleID == sampleIDs[[a]]) %>%
                    mutate(sampleID = NULL)
            })
            names(cbList) <- sampleIDs
        }
        # Return the SimpleList
        bcbio <- list(
            cellularBarcodes = cbList) %>%
            as("SimpleList")
    }

    # bcbioSingleCell ====
    new("bcbioSingleCell",
        SummarizedExperiment(
            assays = assays,
            rowData = rowData,
            colData = colData,
            metadata = metadata),
        bcbio = bcbio)
}



# Methods ====
#' @rdname subset
#' @export
setMethod(
    "[",
    signature(x = "bcbioSingleCell",
              i = "ANY",
              j = "ANY"),
    function(x, i, j, ..., drop = FALSE) {
        .subset(x, i, j, ..., drop)
    })
