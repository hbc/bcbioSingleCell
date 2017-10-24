#' Sanitize Markers Output
#'
#' @rdname sanitizeMarkers
#' @name sanitizeMarkers
#'
#' @inheritParams AllGenerics
#'
#' @param markers [Seurat::FindAllMarkers()] return.
#'
#' @return [grouped_df], grouped by cluster ident.
NULL



# Constructors ====
#' @importFrom basejump camel
#' @importFrom dplyr arrange everything group_by left_join mutate rename select
#' @importFrom rlang !! sym
#' @importFrom tibble as_tibble remove_rownames
.sanitizeMarkersSeurat <- function(
    object,
    markers) {
    # Check for original `Seurat::FindAllMarkers()` return
    seuratCols <- c("p_val", "avg_diff", "pct.1", "pct.2", "cluster", "gene")
    if (all(seuratCols %in% colnames(markers))) {
        message("Sanitizing Seurat markers data.frame")
        # Rename specific columns prior to camelCase
        if ("p_val" %in% colnames(markers)) {
            markers <- rename(markers, pvalue = .data[["p_val"]])
        }
        if ("gene" %in% colnames(markers)) {
            markers <- rename(markers, symbol = .data[["gene"]])
        }
    }

    # Sanitize column names into lowerCamelCase
    markers <- camel(markers, strict = FALSE)
    # Check for ensgene and add from `gene2symbol` if necessary
    if (!"ensgene" %in% colnames(markers)) {
        message("Adding missing Ensembl gene identifiers")
        # Use this data frame to convert the gene symbols ("gene" column) back
        # to Ensembl gene identifiers ("ensgene" column)
        gene2symbol <- data.frame(
            ensgene = names(rownames(slot(object, "raw.data"))),
            symbol = rownames(slot(object, "raw.data"))
        )
        markers <- left_join(markers, gene2symbol, by = "symbol")
        # Ensure all the symbols matched
        if (any(is.na(markers[["ensgene"]]))) {
            stop("gene2symbol failed to match all genes", call. = FALSE)
        }
    }
    # Check for annotable annotations and add if necessary
    if (!"biotype" %in% colnames(markers)) {
        message("Joining annotable metadata")
        # Drop the symbols from annotable before join to avoid mismatch
        annotable <- slot(object, "misc") %>%
            .[["bcbio"]] %>%
            .[["annotable"]] %>%
            mutate(symbol = NULL)
        markers <- left_join(markers, annotable, by = "ensgene")
    }

    # Ensure that required columns are present
    requiredCols <- c(
        "biotype",
        "cluster",
        "description",
        "ensgene",
        "pvalue",
        "symbol")
    if (!all(requiredCols %in% colnames(markers))) {
        stop(paste(
            "Marker data.frame must contain:",
            toString(requiredCols)
        ), call. = FALSE)
    }

    markers %>%
        remove_rownames() %>%
        as_tibble() %>%
        camel(strict = FALSE) %>%
        select(c("cluster", "symbol"), everything()) %>%
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
