#' @name extract
#' @inherit base::Extract title params references
#' @author Michael Steinbaugh
#'
#' @description Extract genes by row and cells by column.
#'
#' @details
#' Refer to [cell2sample()] and [selectSamples()] if sample-level extraction
#' is desired. Note that `sampleID` is slotted into [colData()] and defines the
#' cell-to-sample mappings.
#'
#' Unfiltered cellular barcode distributions for the entire dataset, including
#' cells not kept in the matrix will be dropped in favor of the `nCount` column
#' of `colData()`.
#'
#' @inheritParams general
#'
#' @return `SingleCellExperiment`.
#'
#' @examples
#' data(indrops_small)
#'
#' cells <- head(colnames(indrops_small), 100L)
#' head(cells)
#' genes <- head(rownames(indrops_small), 100L)
#' head(genes)
#'
#' ## Subset by cell identifiers.
#' indrops_small[, cells]
#'
#' ## Subset by genes.
#' indrops_small[genes, ]
#'
#' ## Subset by both genes and cells.
#' indrops_small[genes, cells]
NULL



.extract.SCE <-  # nolint
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
        metadata[["subset"]] <- TRUE

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

        # Return ---------------------------------------------------------------
        new(
            Class = class(x)[[1L]],
            makeSingleCellExperiment(
                assays = assays(sce),
                rowRanges <- rowRanges(sce),
                colData <- colData,
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
    definition = .extract.SCE
)



#' @rdname extract
#' @export
setMethod(
    "[",
    signature(
        x = "CellRanger",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    definition = .extract.SCE
)
