#' Show an object
#'
#' @name show
#' @author Michael Steinbaugh
#' @note Updated 2022-05-09.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return Console output.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioSingleCell ====
#' show(bcb)
NULL



## Updated 2019-07-24.
.showHeader <- function(object, version = NULL) {
    cat(paste(class(object), version), sep = "\n")
}



## Using the same internal method for bcbioSingleCell and CellRanger.
## Updated 2019-08-08.
`show,bcbioSingleCell` <- # nolint
    function(object) {
        validObject(object)
        ## Metadata.
        m <- metadata(object)
        ## Row ranges metadata.
        rrm <- metadata(rowRanges(object))
        .showHeader(object, version = m[["version"]])
        filtered <- isSubset("filterCells", names(m))
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
            filtered = filtered
        ))
        ## Extend the SingleCellExperiment method.
        sce <- as(object, "SingleCellExperiment")
        cat(capture.output(show(sce)), sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature(object = "bcbioSingleCell"),
    definition = `show,bcbioSingleCell`
)
