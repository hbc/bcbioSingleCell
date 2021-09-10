#' Extract or replace parts of an object
#'
#' Extract genes by row and cells by column.
#'
#' Refer to `cell2sample()` and `selectSamples()` if sample-level extraction is
#' desired. Note that `sampleId` is slotted into `colData` and defines the
#' cell-to-sample mappings.
#'
#' Unfiltered cellular barcode distributions for the entire dataset, including
#' cells not kept in the matrix will be dropped in favor of the `nCount` column
#' of `colData()`.
#'
#' @name extract
#' @author Michael Steinbaugh
#' @inherit base::Extract params references
#' @note Updated 2021-09-10.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `bcbioSingleCell`.
#'
#' @examples
#' ## bcbioSingleCell ====
#' data(bcb)
#'
#' cells <- head(colnames(bcb))
#' head(cells)
#' genes <- head(rownames(bcb))
#' head(genes)
#'
#' ## Subset by cell identifiers.
#' bcb[, cells]
#'
#' ## Subset by genes.
#' bcb[genes, ]
#'
#' ## Subset by both genes and cells.
#' bcb[genes, cells]
NULL



## Updated 2019-08-20.
`extract,bcbioSingleCell` <-  # nolint
    function(x, i, j, ..., drop = FALSE) {
        validObject(x)
        assert(identical(drop, FALSE))
        ## Genes (rows).
        if (missing(i)) {
            i <- seq_len(nrow(x))
        }
        ## Cells (columns).
        if (missing(j)) {
            j <- seq_len(ncol(x))
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
            return(x)
        }
        ## Metadata ------------------------------------------------------------
        metadata <- metadata(sce)
        if (isTRUE(subset)) {
            metadata[["cellularBarcodes"]] <- NULL
            metadata[["filterCells"]] <- NULL
            metadata[["filterGenes"]] <- NULL
            metadata[["subset"]] <- TRUE
        }
        metadata <- Filter(f = Negate(is.null), x = metadata)
        metadata(sce) <- metadata
        ## Return --------------------------------------------------------------
        sce <- droplevels(sce)
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
