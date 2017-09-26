.validMarkers <- function(object) {
    if (!is(object, "grouped_df")) {
        stop("Object must be 'grouped_df' class")
    }
    if (attributes(object)[["vars"]] != "cluster") {
        stop("Object must be grouped by 'cluster'")
    }
    if (!"ensgene" %in% colnames(object)) {
        stop("Object missing 'ensgene' column")
    }
    TRUE
}
