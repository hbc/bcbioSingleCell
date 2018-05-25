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
#' @inherit bcbioBase::sampleData
#'
#' @importFrom bcbioBase sampleData sampleData<-
#'
#' @inheritParams general
#'
#' @examples
#' # SingleCellExperiment ====
#' x <- cellranger_small
#' sampleData(x, clean = FALSE) %>% glimpse()
#' sampleData(x, clean = TRUE) %>% glimpse()
#'
#' # Assignment support
#' sampleData(x)[["batch"]] <- 1L
#' sampleData(x) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("SingleCellExperiment"),
    function(
        object,
        clean = TRUE,
        interestingGroups,
        return = c("DataFrame", "data.frame", "kable")
    ) {
        assert_is_a_bool(clean)
        return <- match.arg(return)

        data <- metadata(object)[["sampleData"]]
        if (is.null(data)) {
            return(NULL)
        }
        assert_is_data.frame(data)

        # Only return factor columns, if desired
        if (isTRUE(clean)) {
            data <- data[, vapply(data, is.factor, logical(1L)), drop = FALSE]
            # Drop remaining blacklisted columns
            setdiff <- setdiff(colnames(data), bcbioBase::metadataBlacklist)
            data <- data[, setdiff, drop = FALSE]
        } else {
            # Include `interestingGroups` column, if not NULL
            if (missing(interestingGroups)) {
                interestingGroups <- bcbioBase::interestingGroups(object)
            }
            if (length(interestingGroups)) {
                data <- uniteInterestingGroups(data, interestingGroups)
            }
        }

        # Arrange rows by `sampleName` column, if defined
        if ("sampleName" %in% colnames(data)) {
            data <- data[order(data[["sampleName"]]), , drop = FALSE]
        }

        # Return
        if (return == "kable") {
            kable(as.data.frame(data), row.names = FALSE)
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
        clean = TRUE,
        interestingGroups,
        return = c("DataFrame", "data.frame", "kable")
    ) {
        return <- match.arg(return)

        # Use stashed sampleData, if created with bcbioSingleCell
        data <- metadata(object)[["sampleData"]]
        # Otherwise, generate on the fly
        if (is.null(data)) {
            data <- slot(object, "meta.data")
            assert_is_data.frame(data)
            # Create priority columns from `orig.ident`, if necessary
            if (!"sampleName" %in% colnames(data)) {
                data[["sampleName"]] <- data[["orig.ident"]]
            }
            # Remove columns that map to cells and aren't unique per sample
            blacklist <- paste(
                c(
                    "^cellularBarcode$",
                    "^orig\\.ident$",
                    "^phase$",
                    "^res[.0-9]+$"
                ),
                collapse = "|"
            )
            data <- data %>%
                remove_rownames() %>%
                .[, !grepl(x = colnames(.), pattern = blacklist)] %>%
                mutate_if(is.character, as.factor) %>%
                select_if(is.factor) %>%
                mutate_all(droplevels) %>%
                unique()
            assert_has_no_duplicates(data[["sampleName"]])
            rownames(data) <- makeNames(data[["sampleName"]], unique = TRUE)
            data
        }

        # Only return factor columns, if desired
        if (isTRUE(clean)) {
            data <- data[, vapply(data, is.factor, logical(1L)), drop = FALSE]
            # Drop remaining blacklisted columns
            setdiff <- setdiff(colnames(data), bcbioBase::metadataBlacklist)
            data <- data[, setdiff, drop = FALSE]
        } else {
            # Include `interestingGroups` column, if not NULL
            if (missing(interestingGroups)) {
                interestingGroups <- bcbioBase::interestingGroups(object)
                if (is.null(interestingGroups)) {
                    interestingGroups <- "sampleName"
                }
            }
            if (length(interestingGroups)) {
                data <- uniteInterestingGroups(data, interestingGroups)
            }
        }

        # Arrange rows by `sampleName` column, if defined
        if ("sampleName" %in% colnames(data)) {
            data <- data[order(data[["sampleName"]]), , drop = FALSE]
        }

        # Return
        if (return == "kable") {
            kable(as.data.frame(data), row.names = FALSE)
        } else {
            as(data, return)
        }
    }
)



# Assignment methods ===========================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "SingleCellExperiment",
        value = "DataFrame"
    ),
    function(object, value) {
        value <- sanitizeSampleData(value)
        metadata(object)[["sampleData"]] <- as.data.frame(value)
        object
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "seurat",
        value = "DataFrame"
    ),
    getMethod(
        "sampleData<-",
        signature(
            object = "SingleCellExperiment",
            value = "DataFrame"
        )
    )
)
