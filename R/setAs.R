#' Coerce Object
#'
#' @rdname coerce
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @param from Class for which the coerce method will perform coercion.
#'
#' @seealso `help(topic = "coerce", package = "methods")`.
#'
#' @examples
#' \dontrun{
#' seurat <- as(bcb, "seurat")
#' }
NULL



# Constructors ====
#' Keep Essential Metadata
#'
#' @inheritParams coerce
#'
#' @return Metadata [list].
#' @noRd
.fromMetadata <- function(from) {
    metadata(from) %>%
        .[c("version",
            "sampleMetadata",
            "filterParams",
            "annotable",
            "gene2symbol",
            "genomeBuild",
            "interestingGroups",
            "organism",
            "uploadDir")]
}



# Methods ====
#' @rdname coerce
#' @name upgrade-bcbioSingleCell
#' @section Upgrade [bcbioSingleCell] to current version:
#' This method adds support for upgrading `bcbioSCDataSet` objects to the latest
#' [bcbioSingleCell] class version. This should be backwards compatible to
#' [bcbioSingleCell] version 0.0.17. Previous objects saved using
#' `bcbioSinglecell` (note case) will likely fail to load with newer versions of
#' the package.
setAs("bcbioSingleCellANY", "bcbioSingleCell", function(from) {
    # Can extract the information from the slotted SingleCellExperiment
    # Add a warning for `bcbioSCFiltered`, since these only contain a subset
    stop("Draft function")
})



#' @rdname coerce
#' @name coerce-bcbioSingleCell-seurat
#' @section [bcbioSingleCell] to [seurat]:
#' Interally, this begins by calling [Seurat::CreateSeuratObject()] without any
#' additional filtering cutoffs, since we already applied them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new [seurat] class object, using [as()] object
#' coercion. Next, global-scaling normalization is applied to the raw counts
#' with [Seurat::NormalizeData()], which (1) normalizes the gene expression
#' measurements for each cell by the total expression, (2) multiplies this by a
#' scale factor (10,000 by default), and (3) log-transforms the result.
#' [Seurat::FindVariableGenes()] is then called, which calculates the average
#' expression and dispersion for each gene, places these genes into bins, and
#' then calculates a z-score for dispersion within each bin. This helps control
#' for the relationship between variability and average expression. Finally, the
#' genes are scaled and centered using the [Seurat::ScaleData()] function.
setAs("bcbioSingleCellANY", "seurat", function(from) {
    # Check for required `filterCells` metadata
    if (is.null(metadata(from)[["filterCells"]])) {
        stop(paste(
            "'filterCells()' must be performed on 'from' object",
            "prior to 'seurat' coercion"
        ))
    }

    counts <- counts(from, gene2symbol = TRUE)
    seurat <- CreateSeuratObject(
        raw.data = counts,
        min.cells = 0L,
        min.genes = 0L,
        meta.data = metrics(from)
    )

    # Integrity checks
    if (!identical(dim(counts), dim(seurat@data))) {
        stop(paste(
            "Unexpected dimension mismatch between",
            "'bcbioSingleCell' and 'seurat' objects"
        ))
    }

    # Complete the initalization steps
    seurat <- seurat %>%
        NormalizeData() %>%
        FindVariableGenes(do.plot = FALSE) %>%
        ScaleData()

    # Stash useful bcbio run metadata into `misc` slot
    seurat@misc[["bcbio"]] <- .fromMetadata(from)

    seurat
})
