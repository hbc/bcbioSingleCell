#' Filter Cells
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filterCells
#' @name filterCells
#' @author Michael Steinbaugh
#'
#' @param minUMIs Minimum number of UMI disambiguated counts per cell.
#' @param minGenes Minimum number of genes detected.
#' @param maxGenes Maximum number of genes detected.
#' @param maxMitoRatio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param minNovelty Minimum novelty score.
#' @param showReport Show summary statistics report and plots.
#'
#' @return [bcbioSCSubset].
#' @export
NULL



# Constructors ====
.filterCells <- function(
    object,
    minUMIs = 1000L,
    minGenes = 500L,
    maxGenes = NULL,
    maxMitoRatio = 0.1,
    minNovelty = 0.8,
    showReport = TRUE) {
    sparseCounts <- assay(object)

    # Cellular barcode count
    ncol(sparseCounts) %>%
        paste("cellular barcodes in dataset") %>%
        message

    # Subset metrics
    metrics <- metrics(object)
    if (!is.null(minUMIs)) {
        metrics <- metrics %>%
            .[.[["nUMI"]] >= minUMIs, ]
    }
    if (!is.null(minGenes)) {
        metrics <- metrics %>%
            .[.[["nGene"]] >= minGenes, ]
    }
    if (!is.null(maxGenes)) {
        metrics <- metrics %>%
            .[.[["nGene"]] <= maxGenes, ]
    }
    if (!is.null(maxMitoRatio)) {
        metrics <- metrics %>%
            .[.[["mitoRatio"]] <= maxMitoRatio, ]
    }
    if (!is.null(minNovelty)) {
        metrics <- metrics %>%
            .[.[["log10GenesPerUMI"]] >= minNovelty, ]
    }
    if (!nrow(metrics)) {
        stop("No cellular barcodes passed filtering")
    }
    message(paste(nrow(metrics), "cellular barcodes passed filtering"))

    # Filter the sparse counts matrix with metrics
    sparseCounts <- sparseCounts[, rownames(metrics)]

    # colData ====
    colData <- colData(object) %>%
        .[rownames(metrics), ]
    rm(metrics)

    # rowData ====
    # Add the rownames back
    rowData <- rowData(object) %>%
        set_rownames(rownames(object))

    # Metadata ====
    oldMeta <- metadata(object) %>%
        .[c("version",
            "pipeline",
            "uploadDir",
            "sampleMetadataFile",
            "sampleMetadata",
            "interestingGroups",
            "genomeBuild",
            "annotable",
            "ensemblVersion",
            "umiType",
            "multiplexedFASTQ",
            "runDate",
            "cbCutoff")]
    newMeta <- list(
        filteringCriteria = c(
            minUMIs = minUMIs,
            minGenes = minGenes,
            maxGenes = maxGenes,
            maxMitoRatio = maxMitoRatio,
            minNovelty = minNovelty))
    metadata <- c(oldMeta, newMeta) %>%
        as("SimpleList")

    # SummarizedExperiment ====
    se <- packageSE(
        sparseCounts,
        colData = colData,
        rowData = rowData,
        metadata = metadata)
    object <- new("bcbioSCSubset", se)

    # Show summary statistics report and plots, if desired
    if (isTRUE(showReport)) {
        mdHeader("Filtering parameters", level = 2L)
        mdList(c(
            paste0("`>= ", minUMIs, "` UMI counts per cell"),
            paste0("`>= ", minGenes, "` genes per cell"),
            paste0("`<= ", maxGenes, "` genes per cell"),
            paste0("`<= ", maxMitoRatio, "` relative mitochondrial abundance"),
            paste0("`>= ", minNovelty, "` novelty score")))

        mdHeader("Filtered metrics plots {.tabset}", level = 2L)

        mdHeader("Reads per cell", level = 3L)
        show(plotReadsPerCell(object))

        mdHeader("Cell counts", level = 3L)
        show(plotCellCounts(object))

        mdHeader("UMI counts per cell", level = 3L)
        show(plotUMIsPerCell(object))

        mdHeader("Genes detected", level = 3L)
        show(plotGenesPerCell(object))

        mdHeader("Mitochondrial counts ratio", level = 3L)
        show(plotMitoRatio(object))

        mdHeader("Novelty", level = 3L)
        show(plotNovelty(object))
    }

    object
}



# Methods ====
#' @rdname filterCells
#' @export
setMethod("filterCells", "bcbioSCDataSet", .filterCells)
