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
setMethod("sampleMetadata", "bcbioSingleCellANY", function(
    object,
    aggregateReplicates = TRUE) {
    meta <- metadata(object)[["sampleMetadata"]] %>%
        as.data.frame()
    if (isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(meta)) {
        meta <- meta %>%
            mutate(sampleName = .data[["sampleNameAggregate"]],
                   sampleID = make.names(.data[["sampleName"]]),
                   sampleNameAggregate = NULL) %>%
            # Here we're keeping the priority and interesting group columns
            # only, so we can collapse down to distinct per sample rows
            dplyr::select(c(
                metaPriorityCols,
                interestingGroups(object)
            )) %>%
            distinct() %>%
            set_rownames(.[["sampleID"]])
    }
    meta
})



#' @rdname sampleMetadata
#' @export
setMethod("sampleMetadata", "seurat", function(object) {
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
