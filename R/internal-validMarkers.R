.validMarkers <- function(object) {
    if (!is(object, "grouped_df")) {
        stop("Object must be 'grouped_df' class", call. = FALSE)
    }
    if (attributes(object)[["vars"]] != "cluster") {
        stop("Object must be grouped by 'cluster'", call. = FALSE)
    }
    if (!"ensgene" %in% colnames(object)) {
        stop("Object missing 'ensgene' column", call. = FALSE)
    }
    TRUE
}
