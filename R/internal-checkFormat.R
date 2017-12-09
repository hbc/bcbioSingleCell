.checkFormat <- function(format) {
    validFormat <- c("ensgene", "symbol")
    if (!format %in% validFormat) {
        stop(paste(
            "'format' must contain:", toString(validFormat)
        ), call. = FALSE)
    }
}
