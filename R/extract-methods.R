#' Extract or Replace Parts of an Object
#'
#' Extract genes by row and cells by column from a `bcbioSingleCell` object.
#'
#' @note Unfiltered cellular barcode distributions for the entire dataset,
#'   including cells not kept in the matrix will be dropped in favor of the
#'   `nCount` column of `colData()`.
#'
#' @name extract
#' @family S4 Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams base::`[`
#' @inheritParams general
#'
#' @seealso
#' - `help("[", "base")`.
#' - `selectSamples()` for subsetting based on sample metadata.
#'
#' @return `bcbioSingleCell`.
#'
#' @examples
#' cells <- head(colnames(indrops_small), 100L)
#' head(cells)
#' genes <- head(rownames(indrops_small), 100L)
#' head(genes)
#'
#' # Subset by cell identifiers.
#' indrops_small[, cells]
#'
#' # Subset by genes.
#' indrops_small[genes, ]
#'
#' # Subset by both genes and cells.
#' indrops_small[genes, cells]
NULL



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
    function(x, i, j, ..., drop = FALSE) {
        validObject(x)

        # Genes
        if (missing(i)) {
            i <- 1L:nrow(x)
        }
        # Cells
        if (missing(j)) {
            j <- 1L:ncol(x)
        }

        if (identical(
            x = c(length(i), length(j)),
            y = dim(x)
        )) {
            message("Returning object unmodified.")
            return(x)
        }

        # Subset using SCE method.
        sce <- as(x, "SingleCellExperiment")
        sce <- sce[i, j, drop = drop]

        genes <- rownames(sce)
        cells <- colnames(sce)

        # Column data ----------------------------------------------------------
        # Ensure factors get releveled.
        colData <- colData(sce) %>%
            as("tbl_df") %>%
            mutate_if(is.character, as.factor) %>%
            mutate_if(is.factor, droplevels) %>%
            as("DataFrame")

        # Metadata -------------------------------------------------------------
        metadata <- metadata(sce)
        metadata <- .updateMetadata(metadata)

        # Drop unfiltered cellular barcode list.
        metadata[["cellularBarcodes"]] <- NULL

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

        metadata[["extract"]] <- TRUE

        # Return ---------------------------------------------------------------
        .new.bcbioSingleCell(
            assays = assays(sce),
            rowRanges <- rowRanges(sce),
            colData <- colData,
            metadata = metadata,
            spikeNames = rownames(sce)[isSpike(sce)]
        )
    }
)
