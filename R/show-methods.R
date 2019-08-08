#' Show an object
#' @name show
#' @author Michael Steinbuagh
#' @inherit methods::show
#' @note Updated 2019-08-08.
#' @examples
#' data(indrops)
#' show(indrops)
NULL



## Updated 2019-07-24.
.showHeader <- function(object, version = NULL) {
    cat(paste(class(object), version), sep = "\n")
}



## Using the same internal method for bcbioSingleCell and CellRanger.
## Updated 2019-08-08.
`show,bcbioSingleCell` <-  # nolint
    function(object) {
        validObject(object)
        ## Metadata.
        m <- metadata(object)
        ## Row ranges metadata.
        rrm <- metadata(rowRanges(object))
        .showHeader(object, version = m[["version"]])
        subset <- "subset" %in% names(m)
        showSlotInfo(list(
            uploadDir = m[["uploadDir"]],
            dates = as.character(c(
                bcbio = m[["runDate"]],
                R = m[["date"]]
            )),
            level = m[["level"]],
            sampleMetadataFile = m[["sampleMetadataFile"]],
            organism = m[["organism"]],
            gffFile = m[["gffFile"]],
            annotationHub = rrm[["annotationHub"]],
            ensemblRelease = rrm[["release"]],
            genomeBuild = rrm[["build"]],
            interestingGroups = m[["interestingGroups"]],
            subset = subset
        ))
        ## Extend the SingleCellExperiment method.
        sce <- as(object, "SingleCellExperiment")
        cat(capture.output(show(sce)), sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("bcbioSingleCell"),
    definition = `show,bcbioSingleCell`
)
