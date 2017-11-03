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
        aggregateReplicates = TRUE) {
        meta <- metadata(object)[["sampleMetadata"]] %>%
            as.data.frame()

        # Aggregate replicates, if necessary
        if (isTRUE(aggregateReplicates) &
            "sampleNameAggregate" %in% colnames(meta)) {
            meta <- meta %>%
                mutate(sampleName = .data[["sampleNameAggregate"]],
                       sampleID = make.names(
                           .data[["sampleName"]], unique = FALSE),
                       sampleNameAggregate = NULL) %>%
                select(unique(c(
                    "sampleID",
                    "sampleName",
                    interestingGroups(object)
                ))) %>%
                distinct()
        } else {
            # Put the priority columns first
            meta <- meta %>%
                select(c("sampleID",
                         "sampleName",
                         "description"),
                       everything())
        }

        interestingGroups <- metadata(object)[["interestingGroups"]]
        meta %>%
            uniteInterestingGroups(interestingGroups) %>%
            # Ensure strings as factors
            mutate_if(is.character, as.factor) %>%
            # Ensure the rownames are set
            set_rownames(.[["sampleID"]])
    })



#' @rdname sampleMetadata
#' @importFrom basejump camel
#' @importFrom dplyr distinct funs mutate_if summarize_all
#' @importFrom magrittr set_rownames
#' @importFrom tibble remove_rownames
#' @export
setMethod(
    "sampleMetadata",
    signature("seurat"),
    function(object) {
        # Check to see if bcbio metadata is stashed in '@misc' slot. This is
        # saved when using our `setAs()` coercion method.
        bcbMeta <- slot(object, "misc") %>%
            .[["bcbio"]] %>%
            .[["sampleMetadata"]]
        if (!is.null(bcbMeta)) {
            message(paste(
                "Using bcbio sample metadata stashed in",
                "'object@misc$bcbio$sampleMetadata'"
            ))
            df <- bcbMeta
        } else {
            message(paste(
                "Attempting to construct sample metadata",
                "from 'seurat@meta.data' slot"
            ))
            df <- slot(object, "meta.data") %>%
                remove_rownames()
            # Drop any columns that appear to be a count (e.g. `nUMI`)
            df <- df[, !grepl(x = colnames(df), pattern = "^n[A-Z]")]
            # Drop any tSNE resolution columns (e.g. `res.0.8`)
            df <- df[, !grepl(x = colnames(df), pattern = "^res\\.")]
            # Drop any blacklisted columns
            blacklist <- c(
                # Seurat
                "orig.ident",
                "percent.mito",
                "G2M.Score",
                "S.Score",
                "Phase",
                # bcbio
                "cellularBarcode",
                "log10GenesPerUMI",
                "mitoRatio")
            df <- df %>%
                .[, !colnames(.) %in% blacklist]

            # Only keep columns containing duplicates. Any remaining metrics
            # columns should be filtered by this step.
            findColsWithDupes <- function(.data) {
                .data %>%
                    as.character() %>%
                    duplicated() %>%
                    any()
            }
            hasDupes <- summarize_all(df, funs(findColsWithDupes))
            df <- df[, as.logical(hasDupes)]

            # Finally, attempt to collapse into distinct rows
            df <- distinct(df) %>%
                camel(strict = FALSE)
        }

        # Ensure strings as factors
        df <- mutate_if(df, is.character, as.factor)
        rownames(df) <- df[["sampleID"]]
        df
    })
