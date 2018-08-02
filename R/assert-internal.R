.isAggregate <- function(object, stop = FALSE) {
    logical <- "aggregate" %in% colnames(object)
    if (
        identical(logical, FALSE) &&
        identical(stop, TRUE)
    ) {
        stop("`aggregate` column is required")
    }
    logical
}



# Check assays for zinbwave weights
.hasZinbwave <- function(object) {
    stopifnot(is(object, "SingleCellExperiment"))
    all(c("normalizedValues", "weights") %in% assayNames(object))
}
