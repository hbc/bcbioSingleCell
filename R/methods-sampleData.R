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
#' @param aggregateReplicates Aggregate technical replicates, if specified. This
#'   function uses values assigned in the `sampleNameAggregate` column of the
#'   internal sample metadata `data.frame`.
#'
#' @return `data.frame`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' sampleData(bcb) %>% glimpse()
#'
#' # seurat ====
#' sampleData(pbmc_small) %>% glimpse()
#' sampleData(seurat) %>% glimpse()
NULL



# Constructors =================================================================
#' Prepare Sample Metadata from Seurat
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
    function(
        object,
        interestingGroups,
        aggregateReplicates = FALSE
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        # Legacy `sampleData` column fix
        if ("sampleData" %in% names(metadata(object))) {
            warn(paste(
                "Legacy `sampleData` slot detected.",
                "Run `updateObject()."
            ))
            metadata(object)[["sampleData"]] <- metadata(object)[["sampleData"]]
        }

        metadata <- metadata(object) %>%
            .[["sampleData"]] %>%
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
                unique()
            if (!identical(nrow(metadata), expected)) {
                abort("Failed to aggregate sample metadata uniquely")
            }
        }
        .sanitizeSampleMetadata(
            metadata = metadata,
            interestingGroups = interestingGroups
        )
    })



#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("seurat"),
    function(
        object,
        interestingGroups
    ) {
        metadata <- bcbio(object, "sampleData")
        if (!is.null(metadata)) {
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
            # Fall back to constructing metadata from cellular barcode info
            if (!.hasSlot(object, "version")) {
                abort("Failed to detect seurat version")
            }
            # Access the metadata
            if (.hasSlot(object, "meta.data")) {
                metadata <- slot(object, "meta.data") %>%
                    .sampleData.seurat()
            } else if (.hasSlot(object, "data.info")) {
                # Legacy support for older seurat objects (e.g. pbmc33k)
                metadata <- slot(object, "data.info") %>%
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
            metadata = metadata,
            interestingGroups = interestingGroups
        )
    })
