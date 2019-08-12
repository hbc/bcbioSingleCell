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
#' @note Updated 2019-08-08.
#'
#' @inheritParams acidroxygen::params
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
        assert(identical(drop, FALSE))

        ## Genes (rows).
        if (missing(i)) {
            i <- 1L:nrow(x)
        }
        ## Cells (columns).
        if (missing(j)) {
            j <- 1L:ncol(x)
        }

        ## Determine whether we should stash subset in metadata.
        if (identical(x = dim(x), y = c(length(i), length(j)))) {
            subset <- FALSE
        } else {
            subset <- TRUE
        }

        ## Subset using SCE method.
        sce <- as(x, "SingleCellExperiment")
        sce <- sce[i, j, drop = drop]

        ## Early return original object, if unmodified.
        if (identical(assay(sce), assay(x))) {
            message("Returning unmodified.")
            return(x)
        }

        ## Row data ------------------------------------------------------------
        ## Ensure factors get releveled, if necessary.
        rowRanges <- relevel(rowRanges(sce))

        ## Column data ---------------------------------------------------------
        ## Ensure factors get releveled, if necessary.
        colData <- relevel(colData(sce))

        ## Metadata ------------------------------------------------------------
        metadata <- metadata(sce)
        metadata[["cellularBarcodes"]] <- NULL
        metadata[["filterCells"]] <- NULL
        metadata[["filterGenes"]] <- NULL
        if (isTRUE(subset)) {
            metadata[["subset"]] <- TRUE
        }
        metadata <- Filter(f = Negate(is.null), x = metadata)

        ## Return --------------------------------------------------------------
        ## Note that this approach will keep `spikeNames` set, if defined.
        rowRanges(sce) <- rowRanges
        colData(sce) <- colData
        metadata(sce) <- metadata
        new(Class = "bcbioSingleCell", sce)
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
