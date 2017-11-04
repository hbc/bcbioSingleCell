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
#' @return [data.frame].
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



# Methods ====
#' @rdname sampleMetadata
#' @importFrom dplyr distinct everything mutate mutate_if select
#' @importFrom magrittr set_rownames
#' @importFrom stringr str_match
#' @export
setMethod(
    "sampleMetadata",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups,
        aggregateReplicates = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        metadata <- metadata(object)[["sampleMetadata"]] %>%
            as.data.frame()
        # Aggregate replicates, if necessary
        if (isTRUE(aggregateReplicates) &
            "sampleNameAggregate" %in% colnames(metadata)) {
            metadata <- metadata %>%
                mutate(sampleName = .data[["sampleNameAggregate"]],
                       sampleID = make.names(
                           .data[["sampleName"]], unique = FALSE),
                       sampleNameAggregate = NULL) %>%
                select(unique(c(metadataPriorityCols, interestingGroups))) %>%
                distinct()
        } else {
            # Put the priority columns first
            metadata <- metadata %>%
                select(c(metadataPriorityCols), everything())
        }
        metadata %>%
            uniteInterestingGroups(interestingGroups) %>%
            # Ensure strings as factors
            mutate_if(is.character, as.factor) %>%
            # Ensure the rownames are set
            set_rownames(.[["sampleID"]])
    })



#' @rdname sampleMetadata
#' @importFrom dplyr mutate_if
#' @importFrom magrittr set_rownames
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
                unique(metadata[["sampleID"]]),
                unique(slot(object, "meta.data")[["sampleID"]])
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
        metadata %>%
            # Ensure strings as factors
            mutate_if(is.character, as.factor) %>%
            uniteInterestingGroups(interestingGroups) %>%
            set_rownames(.[["sampleID"]])
    })
