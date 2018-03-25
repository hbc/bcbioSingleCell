#' Show an Object
#'
#' @name show
#' @author Michael Steinbuagh
#'
#' @inherit methods::show
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' show(bcb_small)
NULL



# Methods ======================================================================
#' @rdname show
#' @export
setMethod(
    "show",
    signature("bcbioSingleCell"),
    function(object) {
        validObject(object)
        cat(
            paste(class(object), metadata(object)[["version"]]),
            paste("samples:", nrow(sampleData(object))),
            paste("cells:", ncol(object)),
            paste0(metadata(object)[["level"]], ": ", nrow(object)),
            paste("organism:", metadata(object)[["organism"]]),
            paste("bcbio dir:", metadata(object)[["uploadDir"]]),
            paste("bcbio date:", metadata(object)[["runDate"]]),
            paste("load date:", metadata(object)[["date"]]),
            sep = "\n"
        )
    }
)
