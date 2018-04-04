#' Column Data
#'
#' @name colData
#' @family S4 Class Definition
#' @author Michael Steinbaugh
#'
#' @importFrom SummarizedExperiment colData
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell ====
#' colData(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' colData(seurat_small) %>% glimpse()
#' colData(pbmc_small) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname colData
#' @export
setMethod(
    "colData",
    signature("bcbioSingleCell"),
    function(x) {
        sampleData <- as(metadata(x)[["sampleData"]], "DataFrame")
        colData <- slot(x, "colData")
        assert_are_disjoint_sets(
            x = colnames(sampleData),
            y = colnames(colData)
        )
        cell2sample <- metadata(x)[["cell2sample"]]
        assert_is_factor(cell2sample)
        sampleID <- DataFrame("sampleID" = cell2sample)
        colData <- cbind(colData, sampleID)
        colData[["rowname"]] <- rownames(colData)
        data <- merge(
            x = colData,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        )
        # Ensure the numeric metrics columns appear first
        data <- data[, unique(c(colnames(colData), colnames(data)))]
        rownames(data) <- data[["rowname"]]
        data[["rowname"]] <- NULL
        # Add `interestingGroups` column
        interestingGroups <- metadata(x)[["interestingGroups"]]
        data <- uniteInterestingGroups(data, interestingGroups)
        data
    }
)



#' @rdname seurat
#' @export
setMethod(
    "colData",
    signature("seurat"),
    function(x) {
        sampleData <- sampleData(x, return = "DataFrame")
        colData <- slot(x, "meta.data")
        assert_is_data.frame(colData)
        colData <- as(colData, "DataFrame")
        assert_are_disjoint_sets(
            x = colnames(sampleData),
            y = colnames(colData)
        )
        colData <- camel(colData)
        cell2sample <- cell2sample(x)
        sampleID <- DataFrame("sampleID" = cell2sample)
        colData <- cbind(colData, sampleID)
        colData[["rowname"]] <- rownames(colData)
        data <- merge(
            x = colData,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        )
        # Ensure the numeric metrics columns appear first
        data <- data[, unique(c(colnames(colData), colnames(data)))]
        rownames(data) <- data[["rowname"]]
        data[["rowname"]] <- NULL
        # Add `interestingGroups` column
        interestingGroups <- interestingGroups(x)
        data <- uniteInterestingGroups(data, interestingGroups)
        data
    }
)
