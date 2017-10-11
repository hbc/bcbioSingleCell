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
            "uploadDir",
            "sampleMetadata",
            "interestingGroups",
            "filterParams",
            "organism",
            "genomeBuild",
            "ensemblVersion",
            "annotable",
            "gene2symbol")]
}



.upgradeFromLegacy <- function(from) {
    # Check for version
    version <- metadata(from)[["version"]]
    if (is.null(version)) {
        stop(paste(
            "Unknown bcbio object version.",
            "Please reload with 'loadSingleCell()'."
        ), call. = FALSE)
    }
    message(paste(
        "Upgrading to",
        packageVersion("bcbioSingleCell"),
        "from",
        version
    ))
    message(paste("Existing metadata:", toString(names(metadata(from)))))

    assays <- assays(from)
    rowData <- rowData(from)
    rownames(rowData) <- slot(from, "NAMES")
    colData <- colData(from)

    metadata <- metadata(from)
    metadata[["originalVersion"]] <- metadata[["version"]]
    metadata[["version"]] <- packageVersion("bcbioSingleCell")
    metadata[["upgradeDate"]] <- Sys.Date()

    # Version-specific modifications ====
    if (version <= package_version("0.0.18")) {
        bcbio <- slot(from, "callers")
        # Remove GTF file, if present (too large)
        metadata[["gtf"]] <- NULL
    } else {
        bcbio <- slot(from, "bcbio")
    }

    se <- SummarizedExperiment(
        assays = assays,
        rowData = rowData,
        colData = colData,
        metadata = metadata)

    # Return updated object ====
    to <- new("bcbioSingleCell", se)
    slot(to, "bcbio") <- bcbio
    # Recalculate the cellular barcode metrics
    colData(to) <- calculateMetrics(to)
    validObject(to)
    to
}



.coerceToSeurat <- function(from) {
    cells <- metadata(from)[["filterCells"]]
    if (is.null(cells)) {
        stop(paste(
            "'filterCells()' must be performed on 'from' object",
            "prior to 'seurat' coercion"
        ), call. = FALSE)
    }

    genes <- metadata(from)[["filterCells"]]

    # Subset the object to only contain filtered cells
    from <- from[, cells]

    # FIXME Pass the filtering values here

    counts <- counts(from, gene2symbol = TRUE)
    seurat <- CreateSeuratObject(
        raw.data = counts,
        min.cells = 0,
        min.genes = 0,
        meta.data = metrics(from, aggregateReplicates = TRUE)
    )

    # Integrity checks
    if (!identical(dim(counts), dim(seurat@raw.data))) {
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
setAs(
    "bcbioSingleCellLegacy",
    "bcbioSingleCell",
    .upgradeFromLegacy)



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
setAs(
    "bcbioSingleCell",
    signature("seurat"),
    .coerceToSeurat)
