#' Sample Metadata
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#'
#' @importFrom basejump sampleMetadata
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
#'     file.path("inst", "extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("inst", "extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' sampleMetadata(bcb) %>% glimpse()
#'
#' # seurat
#' sampleMetadata(seurat) %>% glimpse()
NULL



# Constructors ====
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
#' @importFrom basejump camel
#' @importFrom dplyr distinct mutate_if select_if
#' @importFrom tibble remove_rownames
.prepareSampleMetadataFromSeurat <- function(metadata) {
    # Assign the required metadata columns from `orig.ident`, if necessary
    if (!all(metadataPriorityCols %in% colnames(metadata))) {
        for (i in seq_len(ncol(metadata))) {
            metadata[[metadataPriorityCols[[i]]]] <- metadata[["orig.ident"]]
        }
    }
    blacklist <- paste(c("Phase", "res\\."), collapse = "|")
    metadata %>%
        remove_rownames() %>%
        .[, !grepl(x = colnames(.), pattern = blacklist)] %>%
        mutate_if(is.character, as.factor) %>%
        select_if(is.factor) %>%
        distinct() %>%
        camel(strict = FALSE)
}



#' @importFrom dplyr everything mutate_if
#' @importFrom magrittr set_rownames
.returnSampleMetadata <- function(metadata, interestingGroups) {
    metadata %>%
        select(c(metadataPriorityCols), everything()) %>%
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        uniteInterestingGroups(interestingGroups) %>%
        set_rownames(.[["sampleID"]])
}



# Methods ====
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
            interestingGroups <- basejump::interestingGroups(object)
        }
        metadata <- metadata(object)[["sampleMetadata"]] %>%
            as.data.frame()
        # Aggregate replicates, if necessary
        if (isTRUE(aggregateReplicates)) {
            .checkAggregate(metadata, stop = TRUE)
            # Get the expected number of rows
            expected <- length(unique(metadata[["sampleNameAggregate"]]))
            metadata <- metadata %>%
                mutate(sampleName = .data[["sampleNameAggregate"]],
                       description = .data[["sampleName"]],
                       sampleID = make.names(
                           .data[["sampleName"]], unique = FALSE),
                       sampleNameAggregate = NULL) %>%
                select(unique(c(metadataPriorityCols, interestingGroups))) %>%
                distinct()
            if (!identical(nrow(metadata), expected)) {
                stop("Failed to aggregate sample metadata uniquely",
                     call. = FALSE)
            }
        }
        .returnSampleMetadata(
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
        metadata <- NULL
        # Attempt to use stashed metadata at `object@misc$bcbio`. This will
        # only exist for seurat class objects created from bcbioSingleCell.
        bcbio <- bcbio(object)
        if (!is.null(bcbio)) {
            metadata <- bcbio[["sampleMetadata"]]
            if (!identical(
                unique(as.character(metadata[["sampleID"]])),
                unique(as.character(slot(object, "meta.data")[["sampleID"]]))
            )) {
                warning(paste(
                    "Dimension mismatch with stashed metadata.",
                    "Using Seurat cellular barcode metadata instead"
                ), call. = FALSE)
                metadata <- NULL
            }
            # Define interesting groups
            if (missing(interestingGroups)) {
                interestingGroups <- bcbio[["interestingGroups"]]
            }
        }
        # Fall back to constructing metadata from cellular barcode info
        if (is.null(metadata)) {
            if (!.hasSlot(object, "version")) {
                warning("Failed to detect seurat version", call. = FALSE)
            }
            if (.hasSlot(object, "meta.data")) {
                metadata <- slot(object, "meta.data") %>%
                    .prepareSampleMetadataFromSeurat()
            } else if (.hasSlot(object, "data.info")) {
                # Legacy support for older seurat objects (e.g. pbmc33k)
                metadata <- slot(object, "data.info") %>%
                    .prepareSampleMetadataFromSeurat()
            } else {
                stop("Failed to detect metadata in seurat object",
                     call. = FALSE)
            }
            # Define interesting groups
            if (missing(interestingGroups)) {
                interestingGroups <- "sampleName"
            }
        }
        .returnSampleMetadata(
            metadata = metadata,
            interestingGroups = interestingGroups)
    })
