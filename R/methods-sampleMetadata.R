#' Sample Metadata
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#'
#' @inheritParams AllGenerics
#' @inheritParams metrics
#'
#' @return [data.frame].
NULL



# Methods ====
#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata",
    signature("bcbioSingleCellANY"),
    function(
        object,
        aggregateReplicates = TRUE) {
        meta <- metadata(object)[["sampleMetadata"]] %>%
            as.data.frame()
        # Check for and assign missing description (deprecate in future update)
        if (!"description" %in% colnames(meta)) {
            if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
                # `description` is missing in some older bcbio objects because we
                # used `fileName` and `sampleName` initially to define the minimal
                # sample metadata. Now `description` is used for multiplexed
                # samples in QC plots.
                meta[["description"]] <- str_match(
                    meta$sampleID,
                    pattern = "^(.+)_[ACGT]+$") %>%
                    .[, 2]
            } else {
                # `sampleName` is from `description` for demultiplexed samples
                meta[["description"]] <- meta[["sampleName"]]
            }
        }
        if (isTRUE(aggregateReplicates) &
            "sampleNameAggregate" %in% colnames(meta)) {
            meta <- meta %>%
                mutate(sampleName = .data[["sampleNameAggregate"]],
                       sampleID = make.names(.data[["sampleName"]]),
                       sampleNameAggregate = NULL) %>%
                # Here we're keeping the sampleName and interesting group
                # columns only, so we can collapse down to distinct per sample
                # rows
                dplyr::select(
                    unique(c("sampleID",
                             "sampleName",
                             interestingGroups(object)))
                ) %>%
                distinct()
        } else {
            # Put the priority columns first
            meta <- meta %>%
                dplyr::select(c("sampleID", "sampleName", "description"),
                              everything())
        }
        meta %>%
            # Ensure the rownames are set
            set_rownames(.[["sampleID"]])
    })



#' @rdname sampleMetadata
#' @export
setMethod(
    "sampleMetadata",
    signature("seurat"),
    function(object) {
        # Check to see if bcbio metadata is stashed in '@misc' slot. This is
        # saved when using our `setAs()` coercion method.
        bcbMeta <- object@misc[["bcbio"]][["sampleMetadata"]]
        if (!is.null(bcbMeta)) {
            message("Using bcbio sample metadata stashed in '@misc$bcbio'")
            df <- bcbMeta
        } else {
            message(paste("Attempting to construct sample metadata",
                          "from 'seurat@meta.data' slot"))
            df <- object@meta.data %>%
                remove_rownames()
            # Drop any columns that appear to be a count (e.g. `nUMI`)
            df <- df[, !str_detect(colnames(df), "^n[A-Z]")]
            # Drop any tSNE resolution columns (e.g. `res.0.8`)
            df <- df[, !str_detect(colnames(df), "^res\\.")]
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
                camel(strict = FALSE) %>%
                set_rownames(.[["sampleID"]])
        }
        df
    })
