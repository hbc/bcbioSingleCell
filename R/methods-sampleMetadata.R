#' Sample Metadata
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#'
#' @importFrom bcbioBase sampleMetadata
#'
#' @inheritParams AllGenerics
#' @inheritParams metrics
#'
#' @param aggregateReplicates Aggregate technical replicates, if specified. This
#'   function uses values assigned in the `sampleNameAggregate` column of the
#'   internal sample metadata [data.frame].
#'
#' @return [data.frame].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' sampleMetadata(bcb) %>% glimpse()
#'
#' # seurat
#' sampleMetadata(seurat) %>% glimpse()
NULL



# Constructors =================================================================
#' Prepare Sample Metadata from Seurat
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param metadata Seurat metadata [data.frame] per cellular barcode. Currently
#'   this is slotted in `object@meta.data` but older versions used
#'   `object@data.info`.
#'
#' @importFrom bcbioBase camel
#' @importFrom dplyr distinct mutate_if select_if
#' @importFrom tibble remove_rownames
.sampleMetadata.seurat <- function(metadata) {  # nolint
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
        distinct() %>%
        camel(strict = FALSE)

    # Check for failure to make rows distinct (by `sampleName`)
    if (any(duplicated(metadata[["sampleName"]]))) {
        abort("Failed to make `sampleName` column unique")
    }

    metadata
}



#' @importFrom dplyr everything mutate_if
#' @importFrom magrittr set_rownames
.sanitizeSampleMetadata <- function(metadata, interestingGroups) {
    metadata %>%
        select(c(metadataPriorityCols), everything()) %>%
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        uniteInterestingGroups(interestingGroups) %>%
        set_rownames(.[["sampleID"]])
}



# Methods ======================================================================
#' @rdname sampleMetadata
#' @importFrom dplyr distinct mutate select
#' @importFrom stringr str_match
#' @export
setMethod(
    "sampleMetadata",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups,
        aggregateReplicates = FALSE) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        metadata <- metadata(object) %>%
            .[["sampleMetadata"]] %>%
            as.data.frame()
        # Aggregate replicates, if necessary
        if (isTRUE(aggregateReplicates)) {
            .checkAggregate(metadata, stop = TRUE)
            # Get the expected number of rows
            expected <- length(unique(metadata[["sampleNameAggregate"]]))
            metadata <- metadata %>%
                mutate(
                    sampleName = .data[["sampleNameAggregate"]],
                    description = .data[["sampleName"]],
                    sampleID = make.names(
                        .data[["sampleName"]], unique = FALSE),
                    sampleNameAggregate = NULL
                ) %>%
                select(unique(c(metadataPriorityCols, interestingGroups))) %>%
                distinct()
            if (!identical(nrow(metadata), expected)) {
                abort("Failed to aggregate sample metadata uniquely")
            }
        }
        .sanitizeSampleMetadata(
            metadata = metadata,
            interestingGroups = interestingGroups)
    })



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata",
    signature("seurat"),
    function(
        object,
        interestingGroups) {
        # Attempt to use stashed metadata. This will only exist for seurat
        # objects created with bcbioSingleCell.
        metadata <- bcbio(object, "sampleMetadata")
        if (!is.null(metadata)) {
            inform("Using bcbio stashed metadata")
            if (!identical(
                unique(as.character(metadata[["sampleID"]])),
                unique(as.character(slot(object, "meta.data")[["sampleID"]]))
            )) {
                abort(paste(
                    "`sampleID` mismatch with `seurat@meta.data`",
                    "Using Seurat cellular barcode metadata instead"
                ))
            }
            # Define interesting groups
            if (missing(interestingGroups)) {
                interestingGroups <- bcbioBase::interestingGroups(object)
            }
        } else {
            inform("Generating sample metadata from `meta.data` slot")
            # Fall back to constructing metadata from cellular barcode info
            if (!.hasSlot(object, "version")) {
                abort("Failed to detect seurat version")
            }
            # Access the metadata
            if (.hasSlot(object, "meta.data")) {
                metadata <- slot(object, "meta.data") %>%
                    .sampleMetadata.seurat()
            } else if (.hasSlot(object, "data.info")) {
                # Legacy support for older seurat objects (e.g. pbmc33k)
                metadata <- slot(object, "data.info") %>%
                    .sampleMetadata.seurat()
            } else {
                abort("Failed to locate metadata in seurat object")
            }
            # Define interesting groups
            if (missing(interestingGroups)) {
                interestingGroups <- "sampleName"
            }
        }
        .sanitizeSampleMetadata(
            metadata = metadata,
            interestingGroups = interestingGroups)
    })
