#' Subset Per Sample
#'
#' @rdname subsetPerSample
#' @name subsetPerSample
#'
#' @param dir Output directory where to save the subset data.
NULL



# Methods ====
#' @rdname subsetPerSample
#' @export
setMethod("subsetPerSample", "bcbioSCFiltered", function(
    object,
    dir = "data") {
    sampleIDs <- sampleMetadata(object) %>%
        pull("sampleID")
    pbsapply(seq_along(sampleIDs), function(a) {
        sampleID <- sampleIDs[[a]]
        subset <- selectSamples(
            object,
            sampleID = sampleID)
        assign(sampleID, subset)
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        save(list = sampleID,
             file = file.path(dir, paste0(sampleID, ".rda")))
        sampleID
    }) %>%
        as.character
})
