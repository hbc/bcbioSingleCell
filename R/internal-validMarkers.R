.validMarkers <- function(object) {
    # Check for `sanitizedMarkers()` return
    if (attributes(object)[["vars"]] != "cluster" |
        !"ensgene" %in% colnames(object)) {
        stop("'sanitizeMarkers()' return is required")
    }
}
