#' Combine Multiple `SingleCellExperiment ` Objects
#'
#' @note We're attempting to make this as strict as possible, requiring:
#'
#' - Rows (genes) across objects must be identical.
#' - Organism and genome build must be identical.
#' - rowRanges metadata must be identical.
#' - colData (cellular barcode metrics) must contain the same columns.
#' - sampleData must contain the same columns.
#'
#' @name combine
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics combine
#' @export
#'
#' @inheritParams BiocGenerics::combine
#'
#' @seealso `help("merge.Matrix", "Matrix.utils")`.
#'
#' @return `SingleCellExperiment`.
#'
#' @examples
#' # For this simple working example, let's duplicate our minimal dataset,
#' # modify the cell barcodes, and combine
#'
#' x <- indrops_small
#' validObject(x)
#' sampleData(x)
#' sampleNames(x)
#'
#' # Here we're faking a distinct replicate, just as an example
#' y <- indrops_small
#' # Modify the cellular barcodes prefix
#' colnames(y) <- gsub(
#'     pattern = "^multiplexed_AAAAAAAA",
#'     replacement = "multiplexed_CCCCCCCC",
#'     x = colnames(y)
#' )
#' # Modify the column data
#' colData(y)[["description"]] <- factor("multiplexed-CCCCCCCC")
#' colData(y)[["sampleID"]] <- factor("multiplexed_CCCCCCCC")
#' colData(y)[["sampleName"]] <- factor("rep_2")
#' colData(y)[["sequence"]] <- factor("GGGGGGGG")
#' colData(y)[["revcomp"]] <- factor("CCCCCCCC")
#' # Modify the sample data
#' sampleData(y)[["sampleName"]] <- factor("rep_2")
#' sampleData(y)[["description"]] <- factor("multiplexed-CCCCCCCC")
#' sampleData(y)[["sequence"]] <- factor("GGGGGGGG")
#' sampleData(y)[["revcomp"]] <- factor("CCCCCCCC")
#' rownames(sampleData(y)) <- "multiplexed_CCCCCCCC"
#' validObject(y)
#' sampleData(y)
#' sampleNames(y)
#'
#' # Combine two SingleCellExperiment objects
#' c <- combine(x, y)
NULL



# Methods ======================================================================
#' @rdname combine
#' @export
setMethod(
    "combine",
    signature(
        x = "SingleCellExperiment",
        y = "SingleCellExperiment"
    ),
    function(x, y, ...) {
        # Coerce the objects to SingleCellExperiment
        x <- as(x, "SingleCellExperiment")
        y <- as(y, "SingleCellExperiment")

        xGenome <- metadata(x)[c("organism", "genomeBuild")]
        yGenome <- metadata(y)[c("organism", "genomeBuild")]
        assert_are_identical(xGenome, yGenome)

        # Currently we're being strict and requiring that the rows (features)
        # are identical, otherwise zero counts may be misleading
        assert_are_identical(rownames(x), rownames(y))

        # Require that there are no duplicate cells
        assert_are_disjoint_sets(colnames(x), colnames(y))

        # Counts ===============================================================
        xCounts <- counts(x)
        yCounts <- counts(y)
        assert_are_identical(class(xCounts), class(yCounts))
        counts <- cbind(xCounts, yCounts)

        # Row data =============================================================
        # Require that the gene annotations are identical (strict)
        xRanges <- rowRanges(x)
        yRanges <- rowRanges(y)
        assert_are_identical(xRanges, yRanges)
        rowRanges <- xRanges

        # Column data ==========================================================
        xColData <- colData(x)
        yColData <- colData(y)
        assert_are_identical(colnames(xColData), colnames(yColData))
        colData <- rbind(xColData, yColData)

        # Sample data ==========================================================
        xSampleData <- metadata(x)[["sampleData"]]
        ySampleData <- metadata(y)[["sampleData"]]
        assert_are_identical(colnames(xSampleData), colnames(ySampleData))
        sampleData <- rbind(xSampleData, ySampleData)

        # cell2sample ==========================================================
        cell2sample <- mapCellsToSamples(
            cells = colnames(counts),
            samples = rownames(sampleData)
        )

        # Metadata (minimal) ===================================================
        m <- metadata(x)
        metadata <- list(
            interestingGroups = m[["interestingGroups"]],
            organism = m[["organism"]],
            genomeBuild = m[["genomeBuild"]],
            ensemblRelease = m[["ensemblRelease"]],
            rowRangesMetadata = m[["rowRangesMetadata"]],
            sampleData = sampleData,
            cell2sample = cell2sample,
            umiType = m[["umiType"]]
        )

        # Return SingleCellExperiment ==========================================
        .new.SingleCellExperiment(
            assays = list(counts = counts),
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    }
)
