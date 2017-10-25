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
    sanitizedMarkers <- .checkSanitizedMarkers(markers, package = "Seurat")
    # Message and return unmodified, if already sanitized
    if (isTRUE(sanitizedMarkers)) {
        message("Markers are already sanitized")
        return(markers)
    }

    # Update legacy columns
    if ("avg_diff" %in% colnames(markers)) {
        message(paste(
            "Renaming legacy 'avg_diff' to 'avg_logFC'",
            "(changed in Seurat v2.1)"
        ))
    }

    # Rename specific columns prior to camelCase
    if ("gene" %in% colnames(markers)) {
        message("Renaming 'gene' to 'symbol'")
        markers <- rename(markers, symbol = .data[["gene"]])
    }
    if ("p_val" %in% colnames(markers)) {
        message(paste(
            "Renaming 'p_val' to 'pvalue'",
            "(matching DESeq2)"
        ))
        markers <- rename(markers, pvalue = .data[["p_val"]])
    }
    if ("p_val_adj" %in% colnames(markers)) {
        message(paste(
            "Renaming 'p_val_adj' to 'padj'",
            "(matching DESeq2)"
        ))
        markers <- rename(markers, padj = .data[["p_val_adj"]])
    }

    # Sanitize column names into lowerCamelCase
    message("Converting columns into camelCase")
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
    if (!"description" %in% colnames(markers)) {
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
        "avgLogFC",     # Seurat v2.1
        "biotype",      # Ensembl annotations
        "cluster",      # Unmodified
        "description",  # Ensembl annotations
        "ensgene",      # Ensembl annotations
        "pvalue",       # Renamed from `p_val`
        "symbol"        # Renamed from `gene`
    )
    if (!all(requiredCols %in% colnames(markers))) {
        stop(paste(
            "Marker data.frame must contain:",
            toString(requiredCols)
        ), call. = FALSE)
    }

    # Return as a tibble
    # Grouped by cluster
    # Arranged by P value (per cluster)
    markers %>%
        remove_rownames() %>%
        as_tibble() %>%
        # Ensure all the annotations added are camelCase
        camel(strict = FALSE) %>%
        select(c("cluster",
                 "ensgene",
                 "symbol",
                 "pct1",
                 "pct2",
                 "avgLogFC",
                 "pvalue"),
               # `padj` should come after this, but isn't in legacy output
               everything()) %>%
        group_by(.data[["cluster"]]) %>%
        # Arrange by P value
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
