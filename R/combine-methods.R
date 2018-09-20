#' Combine Multiple Objects
#'
#' @name combine
#' @family S4 Functions
#' @author Michael Steinbaugh
#' @importFrom BiocGenerics combine
#' @importMethodsFrom basejump combine
#' @export
#'
#' @inherit basejump::combine
#'
#' @seealso
#' - [BiocGenerics::combine()].
#' - `help("merge.Matrix", "Matrix.utils")`.
#'
#' @return `SingleCellExperiment`.
#'
#' @examples
#' # For this simple working example, let's duplicate our minimal dataset,
#' # modify the cell barcodes, and combine.
#'
#' x <- indrops_small
#' validObject(x)
#' sampleData(x)
#' sampleNames(x)
#'
#' # Here we're faking a distinct replicate, just as an example.
#' y <- indrops_small
#' # Modify the cellular barcodes prefix.
#' colnames(y) <- gsub(
#'     pattern = "^multiplexed_AAAAAAAA",
#'     replacement = "multiplexed_CCCCCCCC",
#'     x = colnames(y)
#' )
#' # Modify the column data.
#' colData(y)[["description"]] <- factor("multiplexed-CCCCCCCC")
#' colData(y)[["sampleID"]] <- factor("multiplexed_CCCCCCCC")
#' colData(y)[["sampleName"]] <- factor("rep_2")
#' colData(y)[["sequence"]] <- factor("GGGGGGGG")
#' colData(y)[["revcomp"]] <- factor("CCCCCCCC")
#' # Modify the sample data.
#' sampleData(y)[["sampleName"]] <- factor("rep_2")
#' sampleData(y)[["description"]] <- factor("multiplexed-CCCCCCCC")
#' sampleData(y)[["sequence"]] <- factor("GGGGGGGG")
#' sampleData(y)[["revcomp"]] <- factor("CCCCCCCC")
#' rownames(sampleData(y)) <- "multiplexed_CCCCCCCC"
#' validObject(y)
#' sampleData(y)
#' sampleNames(y)
#'
#' # Combine two SingleCellExperiment objects.
#' c <- combine(x, y)
#' print(c)
#' sampleNames(c)
NULL



#' @rdname combine
#' @export
setMethod(
    f = "combine",
    signature = signature(
        x = "SingleCellExperiment",
        y = "SingleCellExperiment"
    ),
    definition = function(
        x,
        y,
        metadata = c(
            "version",
            "pipeline",
            "level",
            "interestingGroups",
            "organism",
            "genomeBuild",
            "ensemblRelease",
            "umiType",
            "gffFile",
            "dataVersions",
            "programVersions"
        )
    ) {
        assert_is_character(metadata)

        # Use our SummarizedExperiment method to combine.
        assert_are_identical(class(x), class(y))
        Class <- "RangedSummarizedExperiment"  # nolint
        rse <- combine(
            x = as(object = x, Class = Class),
            y = as(object = y, Class = Class)
        )
        counts <- counts(rse)
        colData <- colData(rse)
        rowRanges <- rowRanges(rse)

        # Update cell2sample mappings ------------------------------------------
        assert_are_set_equal(
            x = colnames(sampleData(x)),
            y = colnames(sampleData(y))
        )
        cols <- intersect(
            x = colnames(sampleData(x)),
            y = colnames(sampleData(y))
        )
        sampleData <- rbind(
            sampleData(x)[, cols, drop = FALSE],
            sampleData(y)[, cols, drop = FALSE]
        )
        cell2sample <- .mapCellsToSamples(
            cells = colnames(counts),
            samples = rownames(sampleData)
        )

        # Metadata -------------------------------------------------------------
        metadata <- metadata(x)[metadata]
        metadata[["cell2sample"]] <- cell2sample

        # Return SingleCellExperiment ------------------------------------------
        .new.SingleCellExperiment(
            assays = list(counts = counts),
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    }
)
