#' Sanitize Markers Output
#'
#' @rdname sanitizeMarkers
#' @name sanitizeMarkers
#'
#' @param object [Seurat::FindAllMarkers()] return.
#'
#' @return [grouped_df], grouped by cluster ident.
NULL



# Constructors ====
.sanitizeMarkersSeurat <- function(
    object,
    gene2symbol) {
    # Check for previously sanitized object and early return, if detected. This
    # function groups by `cluster` ident and ensures Ensembl gene identifiers
    # are added in the `ensgene` column, so this step checks for these.
    if (is(object, "grouped_df") &
        attributes(object)[["vars"]] == "cluster" &
        "ensgene" %in% colnames(object)) {
        return(object)
    }

    # Check for original `Seurat::FindAllMarkers()` return
    seuratCols <- c("p_val", "avg_diff", "pct.1", "pct.2", "cluster", "gene")
    if (all(seuratCols %in% colnames(object))) {
        message("Sanitizing Seurat markers data.frame")
        # Rename specific columns prior to camelCase
        if ("p_val" %in% colnames(object)) {
            object <- dplyr::rename(object, pvalue = .data[["p_val"]])
        }
        if ("gene" %in% colnames(object)) {
            object <- dplyr::rename(object, symbol = .data[["gene"]])
        }
    }

    # Sanitize column names into lowerCamelCase
    object <- camel(object)
    # Check for ensgene and add from `gene2symbol` if necessary
    if (!"ensgene" %in% colnames(object)) {
        message("Adding missing Ensembl gene identifiers")
        # Use this data frame to convert the gene symbols ("gene" column) back to
        # Ensembl gene identifiers ("ensgene" column)
        gene2symbol <- data.frame(
            ensgene = names(rownames(object@raw.data)),
            symbol = rownames(object@raw.data)
        )
        object <- left_join(object, gene2symbol, by = "symbol")
        # Ensure all the symbols matched
        if (any(is.na(object[["ensgene"]]))) {
            stop("gene2symbol failed to match all genes")
        }
    }

    # Ensure that required columns are present
    requiredCols <- c("cluster", "ensgene", "pvalue", "symbol")
    if (!all(requiredCols %in% colnames(object))) {
        stop(paste("Marker data.frame must contain:", toString(requiredCols)))
    }

    object %>%
        remove_rownames() %>%
        as("tibble") %>%
        camel() %>%
        dplyr::select(c("cluster", "symbol"), everything()) %>%
        group_by(.data[["cluster"]]) %>%
        arrange(!!sym("pvalue"), .by_group = TRUE)
}



# Methods ====
#' @rdname sanitizeMarkers
#' @export
setMethod(
    "sanitizeMarkers",
    signature(object = "seurat",
              markers = "data.frame"),
    .sanitizeMarkersSeurat)
