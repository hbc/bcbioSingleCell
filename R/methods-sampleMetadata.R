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
#' @importFrom basejump camel
#' @importFrom dplyr distinct mutate_if select_if
#' @importFrom tibble remove_rownames
.prepareSampleMetadataFromSeurat <- function(object) {
    df <- slot(object, "meta.data") %>%
        as.data.frame() %>%
        remove_rownames()

    # Check for required columns
    if (!all(metadataPriorityCols %in% colnames(df))) {
        stop(paste(
            "Required columns:", toString(metadataPriorityCols)
        ), call. = FALSE)
    }

    # Remove cell cycle regression and resolution columns
    blacklist <- paste(c("Phase", "res\\."), collapse = "|")

    df <- df %>%
        .[, !grepl(x = colnames(.), pattern = blacklist)] %>%
        mutate_if(is.character, as.factor) %>%
        select_if(is.factor) %>%
        distinct() %>%
        camel(strict = FALSE)

    # Check for distinct failure
    if (any(duplicated(df[["sampleID"]]))) {
        stop("Failed to make sample metadata rows distinct", call. = TRUE)
    }

    df
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
            interestingGroups <-
                metadata(object)[["interestingGroups"]]
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
                select(unique(c(
                    "sampleID", "sampleName", interestingGroups
                ))) %>%
                distinct()
        } else {
            # Put the priority columns first
            metadata <- metadata %>%
                select(c("sampleID",
                         "sampleName",
                         "description"),
                       everything())
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
        # Check for stashed metadata in `misc` slot
        bcbio <- bcbio(object)
        if (!is.null(bcbio)) {
            if (missing(interestingGroups)) {
                interestingGroups <- bcbio[["interestingGroups"]]
            }
            df <- bcbio[["sampleMetadata"]]
            if (!identical(
                unique(df[["sampleID"]]),
                unique(slot(object, "meta.data")[["sampleID"]])
            )) {
                stop("Dimension mismatch with stashed metadata",
                     call. = FALSE)
            }
        } else {
            if (missing(interestingGroups)) {
                interestingGroups <- "sampleName"
            }
            df <- .prepareSampleMetadataFromSeurat(object)
        }
        df %>%
            # Ensure strings as factors
            mutate_if(is.character, as.factor) %>%
            uniteInterestingGroups(interestingGroups) %>%
            set_rownames(.[["sampleID"]])
    })
