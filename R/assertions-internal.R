## Updated 2019-07-24.
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
