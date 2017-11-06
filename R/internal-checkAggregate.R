.checkAggregate <- function(object, stop = FALSE) {
    logical <- "sampleNameAggregate" %in% colnames(object)
    if (identical(logical, FALSE) &
        identical(stop, TRUE)) {
        stop("'sampleNameAggregate' column is required",
             call. = FALSE)
    }
    logical
}
