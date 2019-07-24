#' Extract or replace parts of an object
#'
#' Extract genes by row and cells by column.
#'
#' Refer to [`cell2sample()`][basejump::cell2sample] and
#' [selectSamples()][basejump::selectSamples] if sample-level extraction is
#' desired. Note that `sampleID` is slotted into `colData` and defines the
#' cell-to-sample mappings.
#'
#' Unfiltered cellular barcode distributions for the entire dataset, including
#' cells not kept in the matrix will be dropped in favor of the `nCount` column
#' of [`colData()`][SummarizedExperiment::colData].
#'
#' @name extract
#' @author Michael Steinbaugh
#' @inherit base::Extract params references
#' @inheritParams basejump::params
#'
#' @return `bcbioSingleCell`.
#'
#' @examples
#' data(indrops)
#'
#' cells <- head(colnames(indrops), 100L)
#' head(cells)
#' genes <- head(rownames(indrops), 100L)
#' head(genes)
#'
#' ## Subset by cell identifiers.
#' indrops[, cells]
#'
#' ## Subset by genes.
#' indrops[genes, ]
#'
#' ## Subset by both genes and cells.
#' indrops[genes, cells]
NULL



## Updated 2019-07-24.
`extract,bcbioSingleCell` <-  # nolint
    function(x, i, j, ..., drop = FALSE) {
        validObject(x)

        ## Genes
        if (missing(i)) {
            i <- 1L:nrow(x)
        }
        ## Cells
        if (missing(j)) {
            j <- 1L:ncol(x)
        }

        if (identical(
            x = c(length(i), length(j)),
            y = dim(x)
        )) {
            return(x)
        }

        ## Subset using SCE method.
        sce <- as(x, "SingleCellExperiment")
        sce <- sce[i, j, drop = drop]

        genes <- rownames(sce)
        cells <- colnames(sce)

        ## Row data ------------------------------------------------------------
        ## Ensure factors get releveled.
        rowData(sce) <- relevel(rowData(sce))

        ## Column data ---------------------------------------------------------
        ## Ensure factors get releveled.
        colData(sce) <- relevel(colData(sce))

        ## Metadata ------------------------------------------------------------
        metadata <- metadata(sce)
        metadata[["subset"]] <- TRUE

        ## Drop unfiltered cellular barcode list.
        metadata[["cellularBarcodes"]] <- NULL

        ## filterCells
        filterCells <- metadata[["filterCells"]]
        if (!is.null(filterCells)) {
            filterCells <- intersect(filterCells, cells)
            metadata[["filterCells"]] <- filterCells
        }

        ## filterGenes
        filterGenes <- metadata[["filterGenes"]]
        if (!is.null(filterGenes)) {
            filterGenes <- intersect(filterGenes, genes)
            metadata[["filterGenes"]] <- filterGenes
        }

        ## Return --------------------------------------------------------------
        new(
            Class = class(x)[[1L]],
            makeSingleCellExperiment(
                assays = assays(sce),
                rowRanges <- rowRanges(sce),
                colData <- colData(sce),
                metadata = metadata,
                spikeNames = rownames(sce)[isSpike(sce)]
            )
        )
    }



#' @rdname extract
#' @export
setMethod(
    "[",
    signature(
        x = "bcbioSingleCell",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    definition = `extract,bcbioSingleCell`
)
