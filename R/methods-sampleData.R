#' Sample Data
#'
#' Metadata in columns describing the samples, which are defined in the
#' rownames. Similar to [colData()], which for `bcbioSingleCell` and
#' `SingleCellExperiment` objects describes cells in the columns, rather than
#' the samples.
#'
#' @name sampleData
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase sampleData
#'
#' @inheritParams general
#' @inheritParams metrics
#'
#' @return `data.frame`.
#'
#' @examples
#' # bcbioSingleCell ====
#' sampleData(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' sampleData(pbmc_small) %>% glimpse()
#' sampleData(seurat_small) %>% glimpse()
NULL



# Constructors =================================================================
#' Prepare Sample Data from Seurat
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param data Seurat metadata `data.frame` per cellular barcode. Currently
#'   this is slotted in `object@meta.data` but older versions used
#'   `object@data.info`.
#'
#' @importFrom basejump camel
#' @importFrom dplyr mutate_if select_if
#' @importFrom tibble remove_rownames
.sampleData.seurat <- function(data) {  # nolint
    # Assign the required metadata columns from `orig.ident`, if necessary
    if (!all(metadataPriorityCols %in% colnames(data))) {
        missing <- setdiff(metadataPriorityCols, colnames(data))
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

    # Check for failure to make rows distinct (by `sampleName`)
    if (any(duplicated(data[["sampleName"]]))) {
        abort("Failed to make `sampleName` column unique")
    }

    rownames(data) <- data[["sampleID"]]
    data
}



# Methods ======================================================================
#' @rdname sampleData
#' @importFrom basejump sanitizeSampleData
#' @importFrom bcbioBase uniteInterestingGroups
#' @export
setMethod(
    "sampleData",
    signature("bcbioSingleCell"),
    function(object, interestingGroups) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        data <- metadata(object)[["sampleData"]]
        data <- uniteInterestingGroups(data, interestingGroups)
        data <- sanitizeSampleData(data)
        assertHasRownames(data)
        data
    }
)



#' @rdname sampleData
#' @importFrom basejump sanitizeSampleData
#' @importFrom bcbioBase uniteInterestingGroups
#' @export
setMethod(
    "sampleData",
    signature("seurat"),
    function(object, interestingGroups) {
        data <- metadata(object)[["sampleData"]]
        # Define interesting groups
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        if (is.null(data)) {
            # Fall back to constructing metadata from cellular barcode info
            if (!.hasSlot(object, "version")) {
                abort("Failed to detect seurat version")
            }
            # Access the metadata
            if (.hasSlot(object, "meta.data")) {
                data <- slot(object, "meta.data") %>%
                    .sampleData.seurat()
            } else if (.hasSlot(object, "data.info")) {
                # Legacy support for older seurat objects (e.g. pbmc33k)
                data <- slot(object, "data.info") %>%
                    .sampleData.seurat()
            } else {
                abort("Failed to locate metadata in seurat object")
            }
            # Define interesting groups
            if (missing(interestingGroups)) {
                interestingGroups <- "sampleName"
            }
        }
        data <- uniteInterestingGroups(data, interestingGroups)
        data <- sanitizeSampleData(data)
        assertHasRownames(data)
        data
    }
)
