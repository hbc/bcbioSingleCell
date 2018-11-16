#' @name show
#' @author Michael Steinbuagh
#' @inherit methods::show
#' @examples
#' data(indrops)
#' show(indrops)
NULL



#' @importFrom methods show
#' @aliases NULL
#' @export
methods::show



.showHeader <- function(object, version = NULL) {
    cat(paste(class(object), version), sep = "\n")
}



# Using the same internal method for bcbioSingleCell and CellRanger.
show.bcbioSingleCell <- function(object) {
    validObject(object)
    # Metadata.
    m <- metadata(object)
    # Row ranges metadata.
    rrm <- metadata(rowRanges(object))
    .showHeader(object, version = m[["version"]])
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
        filtered = .isFiltered(object)
    ))
    # Extend the SingleCellExperiment method.
    sce <- as(object, "SingleCellExperiment")
    cat(capture.output(show(sce)), sep = "\n")
}



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("bcbioSingleCell"),
    definition = show.bcbioSingleCell
)
