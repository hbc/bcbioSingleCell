#' Filter Cells
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filterCells
#' @name filterCells
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @inheritParams metrics
#' @param minUMIs Minimum number of UMI disambiguated counts per cell.
#' @param minGenes Minimum number of genes detected.
#' @param maxGenes Maximum number of genes detected.
#' @param maxMitoRatio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param minNovelty Minimum novelty score.
#' @param showReport Show summary statistics report and plots.
#' @param headerLevel RMarkdown header level, if `showReport = TRUE`.
#'
#' @note This operation can be re-run on a [bcbioSingleCell] object that has
#'   been previously filtered. This function simply saves `filterCells` and
#'   `filterParams` vectors into the [metadata()] slot.
#'
#' @return [bcbioSingleCell] object, with filtering information slotted into
#'   [metadata()] as `filterCells` and `filterParams`.
#'
#' @examples
#' \dontrun{
#' data(bcb)
#' bcb <- filterCells(bcb)
#' }
NULL



# Constructors ====
.filterCells <- function(
    object,
    aggregateReplicates = FALSE,
    minUMIs = 500,
    minGenes = 500,
    maxGenes = Inf,
    maxMitoRatio = 0.1,
    minNovelty = 0.8,
    showReport = TRUE,
    headerLevel = 2) {
    # Cellular barcode count
    ncol(object) %>%
        paste("cellular barcodes in dataset") %>%
        message()

    # Use the metrics data.frame to identify filtered cells
    metrics <- metrics(
        object,
        filterCells = FALSE,
        aggregateReplicates = aggregateReplicates) %>%
        as.data.frame() %>%
        rownames_to_column("cell")

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

    # Metadata ====
    metadata(object)[["filterCells"]] <- metrics[["cell"]]
    metadata(object)[["filterParams"]] <- c(
        minUMIs = minUMIs,
        minGenes = minGenes,
        maxGenes = maxGenes,
        maxMitoRatio = maxMitoRatio,
        minNovelty = minNovelty
    )

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

        mdHeader(
            "Filtered metrics plots",
            level = headerLevel,
            tabset = TRUE,
            asis = TRUE)

        # Reads per cell currently only supported for bcbio runs
        if (metadata(object)[["pipeline"]] == "bcbio") {
            mdHeader(
                "Reads per cell",
                level = headerLevel + 1,
                asis = TRUE)
            plotReadsPerCell(
                object,
                filterCells = TRUE,
                aggregateReplicates = aggregateReplicates) %>%
                show()
        }

        mdHeader(
            "Cell counts",
            level = headerLevel + 1,
            asis = TRUE)
        plotCellCounts(
            object,
            filterCells = TRUE,
            aggregateReplicates = aggregateReplicates) %>%
            show()

        mdHeader(
            "UMI counts per cell",
            level = headerLevel + 1,
            asis = TRUE)
        plotUMIsPerCell(
            object,
            filterCells = TRUE,
            aggregateReplicates = aggregateReplicates) %>%
            show()

        mdHeader(
            "Genes detected",
            level = headerLevel + 1,
            asis = TRUE)
        plotGenesPerCell(
            object,
            filterCells = TRUE,
            aggregateReplicates = aggregateReplicates) %>%
            show()

        mdHeader(
            "UMIs vs. genes",
            level = headerLevel + 1,
            asis = TRUE)
        plotUMIsVsGenes(
            object,
            filterCells = TRUE,
            aggregateReplicates = aggregateReplicates) %>%
            show()

        mdHeader(
            "Mitochondrial counts ratio",
            level = headerLevel + 1,
            asis = TRUE)
        plotMitoRatio(
            object,
            filterCells = TRUE,
            aggregateReplicates = aggregateReplicates) %>%
            show()

        mdHeader(
            "Novelty",
            level = headerLevel + 1,
            asis = TRUE)
        plotNovelty(
            object,
            filterCells = TRUE,
            aggregateReplicates = aggregateReplicates) %>%
            show()
    }

    object
}



# Methods ====
#' @rdname filterCells
#' @export
setMethod("filterCells", "bcbioSingleCell", .filterCells)



#' @rdname filterCells
#' @export
setMethod("filterCells", "bcbioSCDataSet", function(object) {
    stop("Upgrade 'bcbioSCDataSet' to 'bcbioSingleCell' class object first")
})
