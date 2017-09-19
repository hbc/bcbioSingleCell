#' Sample Metadata
#'
#' @rdname sampleMetadata
#' @name sampleMetadata
#'
#' @inheritParams AllGenerics
#' @param kable Return the output as [kable] instead of [data.frame].
#'
#' @return [data.frame], unless `kable = TRUE`, which will return a [kable].
NULL



# Constructors ====
# This will extract the stashed `sampleMetadata` data.frame
.bcbioSampleMetadata <- function(object, kable = FALSE) {
    object %>%
        metadata %>%
        .[["sampleMetadata"]] %>%
        .returnSampleMetadata(kable = kable)
}



.returnSampleMetadata <- function(object, kable) {
    df <- as.data.frame(object)
    if (isTRUE(kable)) {
        kable(df, caption = "Sample metadata", row.names = FALSE)
    } else {
        df
    }
}



# Methods ====
#' @rdname sampleMetadata
#' @export
setMethod("sampleMetadata", "bcbioSCDataSet", .bcbioSampleMetadata)



#' @rdname sampleMetadata
#' @export
setMethod("sampleMetadata", "bcbioSCFiltered", .bcbioSampleMetadata)



#' @rdname sampleMetadata
#' @export
setMethod("sampleMetadata", "seurat", function(
    object,
    kable = FALSE) {
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
            remove_rownames
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
                as.character %>%
                duplicated %>%
                any
        }
        hasDupes <- summarize_all(df, funs(findColsWithDupes))
        df <- df[, as.logical(hasDupes)]

        # Finally, attempt to collapse into distinct rows
        df <- distinct(df) %>%
            set_rownames(.[["sampleID"]])
    }
    .returnSampleMetadata(df, kable = kable)
})
