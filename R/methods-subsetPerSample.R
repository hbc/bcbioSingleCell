#' Subset Per Sample
#'
#' @rdname subsetPerSample
#' @name subsetPerSample
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @param minCells Minimum number of cells required per sample.
#' @param dir Output directory where to save the subset data.
NULL



# Methods ====
#' @rdname subsetPerSample
#' @export
setMethod("subsetPerSample", "bcbioSCFiltered", function(
    object,
    minCells = 200L,
    dir = "data") {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    sampleIDs <- sampleMetadata(object) %>%
        pull("sampleID")
    pbsapply(seq_along(sampleIDs), function(a) {
        sampleID <- sampleIDs[[a]]
        subset <- selectSamples(
            object,
            sampleID = sampleID)
        # Skip if subset doesn't have enough cells
        if (dim(subset)[[2L]] < minCells) {
            warning(paste(sampleID, "didn't pass minimum cell cutoff"),
                    call. = FALSE)
            return(NA)
        }
        assign(sampleID, subset)
        save(list = sampleID,
             file = file.path(dir, paste0(sampleID, ".rda")))
        sampleID
    }) %>%
        na.omit() %>%
        as.character()
})
