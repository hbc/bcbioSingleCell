#' Sample Data
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
#' @param metadata Seurat metadata `data.frame` per cellular barcode. Currently
#'   this is slotted in `object@meta.data` but older versions used
#'   `object@data.info`.
#'
#' @importFrom basejump camel
#' @importFrom dplyr mutate_if select_if
#' @importFrom tibble remove_rownames
.sampleData.seurat <- function(metadata) {  # nolint
    # Assign the required metadata columns from `orig.ident`, if necessary
    if (!all(metadataPriorityCols %in% colnames(metadata))) {
        missing <- setdiff(metadataPriorityCols, colnames(metadata))
        for (i in seq_along(missing)) {
            metadata[[missing[[i]]]] <- metadata[["orig.ident"]]
        }
    }

    blacklist <- paste(c(
        "cellularBarcode",
        "orig.ident",
        "Phase",
        "^res\\.[0-9]"
    ), collapse = "|")

    metadata <- metadata %>%
        remove_rownames() %>%
        .[, !grepl(x = colnames(.), pattern = blacklist)] %>%
        mutate_if(is.character, as.factor) %>%
        select_if(is.factor) %>%
        unique() %>%
        camel(strict = FALSE)

    # Check for failure to make rows distinct (by `sampleName`)
    if (any(duplicated(metadata[["sampleName"]]))) {
        abort("Failed to make `sampleName` column unique")
    }

    metadata
}



#' @importFrom dplyr everything mutate_all
#' @importFrom magrittr set_rownames
.sanitizeSampleMetadata <- function(metadata, interestingGroups) {
    metadata %>%
        select(c(metadataPriorityCols), everything()) %>%
        mutate_all(as.factor) %>%
        mutate_all(droplevels) %>%
        uniteInterestingGroups(interestingGroups) %>%
        set_rownames(.[["sampleID"]])
}



# Methods ======================================================================
#' @rdname sampleData
#' @importFrom dplyr distinct mutate select
#' @importFrom stringr str_match
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
        .sanitizeSampleMetadata(
            metadata = data,
            interestingGroups = interestingGroups
        )
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("seurat"),
    function(object, interestingGroups) {
        # FIXME Add `metadata()` accessor support to seurat
        data <- bcbio(object, "metadata")[["sampleData"]]
        if (!is.null(data)) {
            # Define interesting groups
            if (missing(interestingGroups)) {
                interestingGroups <- bcbioBase::interestingGroups(object)
            }
        } else {
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
        .sanitizeSampleMetadata(
            metadata = data,
            interestingGroups = interestingGroups
        )
    }
)
