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
        "Upgrading from", version, "to", packageVersion("bcbioSingleCell")
    ))
    message(paste(
        "Existing metadata:", toString(names(metadata(from)))
    ))

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
#' @importFrom Seurat CreateSeuratObject FindVariableGenes NormalizeData
#'   ScaleData
#'
#' @return [seurat] object.
.coerceToSeurat <- function(from) {
    # Require that technical replicates are aggregated
    if (any(grepl(x = metadata(from)[["sampleMetadata"]][["sampleID"]],
                  pattern = lanePattern)) |
        "sampleNameAggregate" %in%
        colnames(metadata(from)[["sampleMetadata"]])) {
        stop(paste("'aggregateReplicates()' required",
                   "to merge technical replicates prior to seurat coercion"
        ), call. = FALSE)
    }

    # Require filtered cells and genes only
    from <- .applyFilterCutoffs(from)
    filterParams <- metadata(from)[["filterParams"]]
    print(filterParams)

    # Create the initial `seurat` object
    rawData <- counts(from, gene2symbol = TRUE)
    metrics <- metrics(from)

    # Ensure that cutoffs are defined
    minGenes <- metadata(from)[["filterParams"]][["minGenes"]]
    if (!is.numeric(minGenes)) {
        stop("'minGenes' is not defined", call. = FALSE)
    }
    message(paste("Minimum genes per cell:", minGenes))
    minCells <- metadata(from)[["filterParams"]][["minCellsPerGene"]]
    # Check for NULL here
    if (!is.numeric(minCells)) {
        stop("'minCells' is not defined", call. = FALSE)
    }
    message(paste("Minimum cells per gene:", minCells))

    # Note here that passing in the `minCells` argument will rescale the number
    # of genes, since we calculated our genes that passed cutoffs based on
    # detection in all cells in the dataset. It makes sense to see the gene
    # count decrease here on a bcbioSingleCell object that has been subset from
    # the main dataset (e.g. using `selectSamples()`).
    seurat <- CreateSeuratObject(
        raw.data = rawData,
        project = "bcbioSingleCell",
        min.cells = minCells,
        min.genes = minGenes,
        # Default for UMI datasets
        is.expr = 0,
        meta.data = metrics) %>%
        NormalizeData(
            object = .,
            normalization.method = "LogNormalize",
            scale.factor = 10000) %>%
        FindVariableGenes(
            object = .,
            mean.function = ExpMean,
            dispersion.function = LogVMR,
            do.plot = FALSE) %>%
        ScaleData(
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
