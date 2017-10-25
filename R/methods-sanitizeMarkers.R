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
    # Check for original `Seurat::FindAllMarkers()` return.
    # These columns are output in an inconsistent format, so we'll sanitize
    # into lowerCamelCase.
    seuratCols <-
        c("avg_diff",   # legacy, now "avg_logFC"
          "avg_logFC",  # renamed in v2.1
          "cluster",
          "gene",       # gene symbol, we'll rename to "symbol"
          "p_val",      # we'll rename to pvalue, matching DESeq2
          "p_val_adj",  # new in v2.1, we'll rename to padj, matching DESeq2
          "pct.1",
          "pct.2")
    if (all(colnames(markers) %in% seuratCols)) {
        message("Original Seurat markers return detected")
        # Rename specific columns prior to camelCase
        if ("gene" %in% colnames(markers)) {
            markers <- rename(markers, symbol = .data[["gene"]])
        }
        if ("p_val" %in% colnames(markers)) {
            markers <- rename(markers, pvalue = .data[["p_val"]])
        }
        if ("p_val_adj" %in% colnames(markers)) {
            markers <- rename(markers, padj = .data[["p_val_adj"]])
        }
    } else {
        stop(paste(
            "Failed to match the original Seurat marker columns.",
            "Has the 'markers' data.frame already been sanitized?",
            "If not, check to see if the Seurat code has changed.",
            "Columns should resemble (legacy included):",
            toString(seuratCols)
        ), call. = FALSE)
    }

    # Sanitize column names into lowerCamelCase
    markers <- camel(markers, strict = FALSE)

    # Check for ensgene and add from `gene2symbol` if necessary
    if (!"ensgene" %in% colnames(markers)) {
        message("Adding stashed Ensembl gene identifiers")
        # When creating the seurat object using the `as()` coercion method, we
        # stashed the gene2symbol identifier mappings both in (1) the rownames
        # of the seurat counts matrix, and (2) the bcbio list of the misc slot.
        # Here we're using the rownames of the counts matrix. The way this
        # works is the gene symbol (`symbol`) is defined as the rowname, and
        # the Ensembl gene identifier (`ensgene`) is defined as the name
        # (`names()`) of the rowname.
        gene2symbol <- data.frame(
            ensgene = names(rownames(slot(object, "raw.data"))),
            symbol = rownames(slot(object, "raw.data"))
        )

        # Check to make sure these values are unique
        if (any(duplicated(gene2symbol[["ensgene"]]))) {
            stop(paste(
                "The seurat object doesn't appear to have",
                "stashed gene identifier mappings"
            ), call. = FALSE)
        }

        markers <- left_join(markers, gene2symbol, by = "symbol")

        # Ensure all the symbols matched
        if (any(is.na(markers[["ensgene"]]))) {
            stop("gene2symbol failed to match all genes", call. = FALSE)
        }
    }

    # Check for annotable annotations and add if necessary
    if (!"biotype" %in% colnames(markers)) {
        message("Adding stashed Ensembl annotations")
        annotable <- slot(object, "misc") %>%
            .[["bcbio"]] %>%
            .[["annotable"]] %>%
            # Drop the symbols from annotable before join to avoid mismatch
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
