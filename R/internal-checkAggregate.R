.checkAggregate <- function(object, stop = FALSE) {
    logical <- "sampleNameAggregate" %in% colnames(object)
    if (identical(logical, FALSE) &
        identical(stop, TRUE)) {
        abort("`sampleNameAggregate` column is required")
    }
    logical
}
