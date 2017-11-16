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
#'
#' @param minUMIs Minimum number of UMI disambiguated counts per cell.
#' @param minGenes Minimum number of genes detected.
#' @param maxGenes Maximum number of genes detected.
#' @param maxMitoRatio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param minNovelty Minimum novelty score.
#' @param minCellsPerGene Include genes with non-zero expression in at least
#'   this many cells.
#' @param showReport Show summary statistics report and plots.
#' @param headerLevel RMarkdown header level, if `showReport = TRUE`.
#' @param destructive Drop low quality cells from object after filtering.
#'
#' @note This operation can be re-run on a [bcbioSingleCell] object that has
#'   been previously filtered. This function simply saves `filterCells` and
#'   `filterParams` vectors into the [metadata()] slot.
#'
#' @seealso [Seurat::CreateSeuratObject()].
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
#' @importFrom basejump mdHeader mdList
#' @importFrom Matrix rowSums
#' @importFrom scales percent
#' @importFrom tibble rownames_to_column
.filterCells <- function(
    object,
    minUMIs = 1000,
    minGenes = 500,
    maxGenes = Inf,
    maxMitoRatio = 0.1,
    minNovelty = 0.8,
    minCellsPerGene = 3,
    showReport = TRUE,
    headerLevel = 2,
    destructive = FALSE) {
    # Ensure that all filter parameters are numeric
    filterParams <- c(
        minUMIs = minUMIs,
        minGenes = minGenes,
        maxGenes = maxGenes,
        maxMitoRatio = maxMitoRatio,
        minNovelty = minNovelty,
        minCellsPerGene = minCellsPerGene
    )
    if (!is.numeric(filterParams)) {
        stop("Filter parameters must all be numeric", call. = FALSE)
    }

    message(paste(ncol(object), "unfiltered cellular barcodes"))

    # Filter low quality cells ====
    # Use the metrics `data.frame` to identify filtered cells
    metrics <- metrics(object, filterCells = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cellID")
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
        stop("No cellular barcodes passed filtering", call. = FALSE)
    }
    filterCells <- metrics[["cellID"]]
    message(paste(
        length(filterCells),
        "filtered cellular barcodes",
        paste0("(", percent(length(filterCells) / ncol(object)), ")")
    ))

    # Filter low expression genes ====
    if (minCellsPerGene > 0) {
        counts <- assay(object)
        numCells <- Matrix::rowSums(counts > 0)
        filterGenes <- names(numCells[which(numCells >= minCellsPerGene)])
    } else {
        filterGenes <- NULL
    }

    # Metadata ====
    metadata(object)[["filterCells"]] <- filterCells
    metadata(object)[["filterGenes"]] <- filterGenes
    metadata(object)[["filterParams"]] <- filterParams

    # Destructive filter ====
    if (isTRUE(destructive)) {
        object <- object[filterGenes, filterCells]
    }

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
            plotReadsPerCell(object, filterCells = TRUE) %>%
                show()
        }

        mdHeader(
            "Cell counts",
            level = headerLevel + 1,
            asis = TRUE)
        plotCellCounts(object, filterCells = TRUE) %>%
            show()

        mdHeader(
            "UMI counts per cell",
            level = headerLevel + 1,
            asis = TRUE)
        plotUMIsPerCell(object, filterCells = TRUE) %>%
            show()

        mdHeader(
            "Genes detected",
            level = headerLevel + 1,
            asis = TRUE)
        plotGenesPerCell(object, filterCells = TRUE) %>%
            show()

        mdHeader(
            "UMIs vs. genes",
            level = headerLevel + 1,
            asis = TRUE)
        plotUMIsVsGenes(object, filterCells = TRUE) %>%
            show()

        mdHeader(
            "Mitochondrial counts ratio",
            level = headerLevel + 1,
            asis = TRUE)
        plotMitoRatio(object, filterCells = TRUE) %>%
            show()

        mdHeader(
            "Novelty",
            level = headerLevel + 1,
            asis = TRUE)
        plotNovelty(object, filterCells = TRUE) %>%
            show()
    }

    object
}



# Methods ====
#' @rdname filterCells
#' @export
setMethod(
    "filterCells",
    signature("bcbioSingleCell"),
    .filterCells)



#' @rdname filterCells
#' @export
setMethod(
    "filterCells",
    signature("bcbioSCDataSet"),
    function(object) {
        stop(paste(
            "Convert 'bcbioSCDataSet' to 'bcbioSingleCell' class.\n",
            "Run this code: bcb <- as(bcb, \"bcbioSingleCell\")"
        ), call. = FALSE)
})
