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
#' Essential Metadata
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom S4Vectors metadata
#'
#' @inheritParams coerce
#'
#' @return Metadata [list].
.essentialMetadata <- function(from) {
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



#' Coerce Legacy bcbio Object to `bcbioSingleCell` class
#'
#' Compatible with old versions created by bcbioSingleCell package.
#' The previous bcbioSinglecell package (note lowercase "c") must be reinstalled
#' to load objects from versions <= 0.0.16.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay assays rowData
#' @importFrom utils packageVersion
#'
#' @return [bcbioSingleCell] object.
.coerceLegacy <- function(from) {
    # Check for version
    version <- metadata(from)[["version"]]
    if (is.null(version)) {
        stop(paste(
            "Unknown bcbio object version.",
            "Please reload with 'loadSingleCell()'."
        ), call. = FALSE)
    }
    message(paste(
        "Upgrading from",
        version,
        "to",
        packageVersion("bcbioSingleCell")
    ))
    message(paste("Existing metadata:", toString(names(metadata(from)))))

    assays <- assays(from)
    rowData <- rowData(from)
    rownames(rowData) <- slot(from, "NAMES")

    metadata <- metadata(from)
    metadata[["originalVersion"]] <- metadata[["version"]]
    metadata[["version"]] <- packageVersion("bcbioSingleCell")
    metadata[["upgradeDate"]] <- Sys.Date()

    # Recalculate the cellular barcode metrics
    colData <- calculateMetrics(
        assay(from),
        metadata[["annotable"]])

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
    validObject(to)
    to
}



#' Coerce `bcbioSingleCell` to `seurat`
#'
#' Last tested against CRAN version 2.0.1
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom BiocGenerics counts
#' @importFrom S4Vectors metadata
#' @importFrom Seurat CreateSeuratObject FindVariableGenes NormalizeData
#'   ScaleData
#'
#' @return [seurat] object.
.coerceToSeurat <- function(from) {
    if (is.null(metadata(from)[["filterParams"]])) {
        stop("'filterCells()' must be run prior to 'seurat' coercion",
             call. = FALSE)
    }

    # Create the initial `seurat` object
    rawData <- counts(from, gene2symbol = TRUE)
    minGenes <- metadata(from)[["filterParams"]][["minGenes"]]
    if (is.null(minGenes)) {
        minGenes <- 0
    }
    minCells <- metadata(from)[["filterParams"]][["minCellsPerGene"]]
    if (is.null(minCells)) {
        minCells <- 0
    }
    metadata <- metrics(
        from,
        aggregateReplicates = FALSE,
        filterCells = FALSE)
    # Add a call to `aggregateReplicates()` here to combine the counts per
    # sample before passing to Seurat, if technical replicates are present?
    seurat <- Seurat::CreateSeuratObject(
        raw.data = rawData,
        project = "bcbioSingleCell",
        min.cells = minCells,
        min.genes = minGenes,
        is.expr = 0,  # Default for UMI datasets
        meta.data = metadata) %>%
        Seurat::NormalizeData(
            object = .,
            normalization.method = "LogNormalize",
            scale.factor = 10000) %>%
        Seurat::FindVariableGenes(
            object = .,
            mean.function = ExpMean,
            dispersion.function = LogVMR,
            do.plot = FALSE) %>%
        Seurat::ScaleData(
            object = .,
            model.use = "linear")

    # Stash useful bcbio run metadata into `misc` slot
    slot(seurat, "misc")[["bcbio"]] <- .essentialMetadata(from)

    print(seurat)
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
    "bcbioSCDataSet",
    signature("bcbioSingleCell"),
    .coerceLegacy)



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
