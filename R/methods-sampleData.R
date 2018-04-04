#' Sample Data
#'
#' Metadata in columns describing the samples, which are defined in the
#' rownames. Similar to [colData()], which for `bcbioSingleCell` and
#' `SingleCellExperiment` objects describes cells in the columns, rather than
#' the samples.
#'
#' @name sampleData
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase sampleData
#'
#' @inheritParams general
#'
#' @return `data.frame`.
#'
#' @examples
#' # bcbioSingleCell ====
#' sampleData(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' sampleData(pbmc_small) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        return = c("DataFrame", "data.frame", "kable")
    ) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        return <- match.arg(return)
        data <- metadata(object)[["sampleData"]]
        assert_is_any_of(data, c("DataFrame", "data.frame"))
        data <- uniteInterestingGroups(data, interestingGroups)
        data <- sanitizeSampleData(data)
        assertHasRownames(data)
        if (return == "kable") {
            blacklist <- c("description", "fileName", "sampleID")
            data %>%
                as.data.frame() %>%
                .[, setdiff(colnames(.), blacklist), drop = FALSE] %>%
                # Ensure `sampleName` is first
                .[, unique(c("sampleName", colnames(.))), drop = FALSE] %>%
                # Arrange by `sampleName`
                .[order(.[["sampleName"]]), , drop = FALSE] %>%
                kable(row.names = FALSE)
        } else {
            as(data, return)
        }
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("seurat"),
    function(
        object,
        interestingGroups,
        return = c("DataFrame", "data.frame")
    ) {
        validObject(object)
        stopifnot(.hasSlot(object, "version"))
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        return <- match.arg(return)
        data <- metadata(object)[["sampleData"]]
        if (is.null(data)) {
            data <- slot(object, "meta.data")
            assert_is_data.frame(data)
            # Create priority columns from `orig.ident`, if necessary
            if (!all(bcbioBase::metadataPriorityCols %in% colnames(data))) {
                missing <- setdiff(
                    x = bcbioBase::metadataPriorityCols,
                    y = colnames(data)
                )
                for (i in seq_along(missing)) {
                    data[[missing[[i]]]] <- data[["orig.ident"]]
                }
            }
            blacklist <- paste(
                c(
                    "cellularBarcode",
                    "orig.ident",
                    "Phase",
                    "^res\\.[0-9]"
                ),
                collapse = "|"
            )
            data <- data %>%
                remove_rownames() %>%
                .[, !grepl(x = colnames(.), pattern = blacklist)] %>%
                mutate_if(is.character, as.factor) %>%
                select_if(is.factor) %>%
                mutate_all(droplevels) %>%
                unique() %>%
                camel()
            assert_has_no_duplicates(data[["sampleName"]])
            rownames(data) <- data[["sampleID"]]
            data
        }
        data <- uniteInterestingGroups(data, interestingGroups)
        data <- sanitizeSampleData(data)
        assertHasRownames(data)
        as(data, return)
    }
)
