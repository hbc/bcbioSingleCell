#' Filter Cells
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filterCells
#' @name filterCells
#'
#' @param minUMIs Minimum number of UMI disambiguated counts per cell.
#' @param minGenes Minimum number of genes detected.
#' @param maxGenes Maximum number of genes detected.
#' @param maxMitoRatio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param minNovelty Minimum novelty score.
#' @param showReport Show summary statistics report and plots.
#' @param headerLevel RMarkdown header level, if `showReport = TRUE`.
#'
#' @return [bcbioSCFiltered].
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
    showReport = TRUE,
    headerLevel = 2L) {
    sparseCounts <- assay(object)

    # Cellular barcode count
    ncol(sparseCounts) %>%
        paste("cellular barcodes in dataset") %>%
        message

    # Subset metrics
    metrics <- metrics(object)
    if (!is.null(minUMIs)) {
        metrics <- metrics %>%
            .[.[["nUMI"]] >= minUMIs, , drop = FALSE]
    }
    if (!is.null(minGenes)) {
        metrics <- metrics %>%
            .[.[["nGene"]] >= minGenes, , drop = FALSE]
    }
    if (!is.null(maxGenes)) {
        metrics <- metrics %>%
            .[.[["nGene"]] <= maxGenes, , drop = FALSE]
    }
    if (!is.null(maxMitoRatio)) {
        metrics <- metrics %>%
            .[.[["mitoRatio"]] <= maxMitoRatio, , drop = FALSE]
    }
    if (!is.null(minNovelty)) {
        metrics <- metrics %>%
            .[.[["log10GenesPerUMI"]] >= minNovelty, , drop = FALSE]
    }
    if (!nrow(metrics)) {
        stop("No cellular barcodes passed filtering")
    }
    message(paste(nrow(metrics), "cellular barcodes passed filtering"))

    # Filter the sparse counts matrix with metrics
    sparseCounts <- sparseCounts[, rownames(metrics), drop = FALSE]

    # colData ====
    colData <- colData(object) %>%
        .[rownames(metrics), , drop = FALSE]
    rm(metrics)

    # rowData ====
    # Add the rownames back
    rowData <- rowData(object) %>%
        set_rownames(rownames(object))

    # Metadata ====
    metadata <- metadata(object)
    metadata[["filterParams"]] <- c(
        minUMIs = minUMIs,
        minGenes = minGenes,
        maxGenes = maxGenes,
        maxMitoRatio = maxMitoRatio,
        minNovelty = minNovelty)

    # SummarizedExperiment ====
    se <- prepareSummarizedExperiment(
        sparseCounts,
        colData = colData,
        rowData = rowData,
        metadata = metadata)
    object <- new("bcbioSCFiltered", se)

    # Show summary statistics report and plots, if desired
    if (isTRUE(showReport)) {
        mdHeader("Filter parameters", level = headerLevel, asis = TRUE)
        mdList(c(
            paste0("`>= ", minUMIs, "` UMI counts per cell"),
            paste0("`>= ", minGenes, "` genes per cell"),
            paste0("`<= ", maxGenes, "` genes per cell"),
            paste0("`<= ", maxMitoRatio, "` relative mitochondrial abundance"),
            paste0("`>= ", minNovelty, "` novelty score")),
            asis = TRUE)

        mdHeader("Filtered metrics plots",
                 level = headerLevel,
                 tabset = TRUE,
                 asis = TRUE)

        # Reads per cell currently only supported for bcbio runs
        if (metadata[["pipeline"]] == "bcbio") {
            mdHeader("Reads per cell", level = headerLevel + 1L, asis = TRUE)
            show(plotReadsPerCell(object))
        }

        mdHeader("Cell counts",
                 level = headerLevel + 1L,
                 asis = TRUE)
        show(plotCellCounts(object))

        mdHeader("UMI counts per cell",
                 level = headerLevel + 1L,
                 asis = TRUE)
        show(plotUMIsPerCell(object))

        mdHeader("Genes detected",
                 level = headerLevel + 1L,
                 asis = TRUE)
        show(plotGenesPerCell(object))

        mdHeader("Mitochondrial counts ratio",
                 level = headerLevel + 1L,
                 asis = TRUE)
        show(plotMitoRatio(object))

        mdHeader("Novelty",
                 level = headerLevel + 1L,
                 asis = TRUE)
        show(plotNovelty(object))
    }

    object
}



# Methods ====
#' @rdname filterCells
#' @export
setMethod("filterCells", "bcbioSCDataSet", .filterCells)
