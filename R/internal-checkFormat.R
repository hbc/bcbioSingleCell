.checkFormat <- function(format) {
    validFormat <- c("ensgene", "symbol")
    if (!format %in% validFormat) {
        abort(paste("`format` must contain:", toString(validFormat)))
    }
}
