#' Sanitize Markers Output
#'
#' @rdname sanitizeMarkers
#' @name sanitizeMarkers
#'
#' @inheritParams general
#'
#' @param markers [Seurat::FindAllMarkers()] return.
#'
#' @return `grouped_df`, grouped by cluster ident.
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#' load(system.file(
#'     "extdata/seuratAllMarkersOriginal.rda",
#'     package = "bcbioSingleCell"
#' ))
#'
#' # seurat ====
#' sanitizeMarkers(
#'     seurat,
#'     markers = seuratAllMarkersOriginal
#' ) %>%
#'     glimpse()
NULL



# Constructors =================================================================
#' @importFrom basejump camel
#' @importFrom dplyr arrange everything group_by left_join mutate select
#' @importFrom rlang !! sym
#' @importFrom tibble as_tibble remove_rownames
.sanitizeMarkersSeurat <- function(
    object,
    markers
) {
    sanitizedMarkers <- .checkSanitizedMarkers(markers, package = "Seurat")
    # Message and return unmodified, if already sanitized
    if (isTRUE(sanitizedMarkers)) {
        inform("Markers are already sanitized")
        return(markers)
    }

    # Update legacy columns
    if ("avg_diff" %in% colnames(markers)) {
        inform(paste(
            "Renaming legacy `avg_diff` column to `avg_logFC`",
            "(changed in Seurat v2.1)"
        ))
    }

    # Rename specific columns prior to camelCase
    if ("gene" %in% colnames(markers)) {
        inform("Renaming `gene` column to `geneName`")
        markers <- dplyr::rename(markers, geneName = .data[["gene"]])
    }
    if ("p_val" %in% colnames(markers)) {
        inform(paste(
            "Renaming `p_val` column to `pvalue`",
            "(matching DESeq2)"
        ))
        markers <- dplyr::rename(markers, pvalue = .data[["p_val"]])
    }
    if ("p_val_adj" %in% colnames(markers)) {
        inform(paste(
            "Renaming `p_val_adj` column to `padj`",
            "(matching DESeq2)"
        ))
        markers <- dplyr::rename(markers, padj = .data[["p_val_adj"]])
    }

    # Sanitize column names into lowerCamelCase
    inform("Converting columns into camelCase")
    markers <- camel(markers)

    # Check for geneID and add from `gene2symbol` if necessary
    if (!"geneID" %in% colnames(markers)) {
        inform("Adding stashed Ensembl gene identifiers")
        # When creating the seurat object using the `as()` coercion method, we
        # stashed the gene2symbol identifier mappings both in (1) the rownames
        # of the seurat counts matrix, and (2) the bcbio list of the misc slot.
        # Here we're using the rownames of the counts matrix. The way this
        # works is the gene name (`geneName`) is defined as the rowname, and
        # the Ensembl gene identifier (`geneID`) is defined as the name
        # (`names()`) of the rowname.
        gene2symbol <- data.frame(
            "geneID" = names(rownames(slot(object, "raw.data"))),
            "geneName" = rownames(slot(object, "raw.data")),
            stringsAsFactors = FALSE
        )

        # Check to make sure these values are unique
        if (any(duplicated(gene2symbol[["geneID"]]))) {
            abort(paste(
                "seurat object doesn't appear to have",
                "stashed gene identifier mappings"
            ))
        }

        markers <- left_join(markers, gene2symbol, by = "geneName")

        # Ensure all the symbols matched
        if (any(is.na(markers[["geneID"]]))) {
            abort("gene2symbol failed to match all genes")
        }
    }

    # Check for rowData annotations and add if necessary
    if (!"description" %in% colnames(markers)) {
        inform("Adding stashed Ensembl annotations")
        rowData <- bcbio(object, "rowData")
        # Drop the symbols from rowData before join to avoid mismatch
        rowData[["geneName"]] <- NULL
        # Ensure Entrez IDs nested as a list get sanitized to string
        if (is.list(rowData[["entrezID"]])) {
            rowData[["entrezID"]] <- vapply(
                rowData[["entrezID"]],
                FUN = toString,
                FUN.VALUE = character(1L)
            )
        }
        markers <- left_join(markers, rowData, by = "geneID")
    }

    # Ensure that required columns are present
    # FIXME This code is duplicated
    # FIXME Switch to ident instead of cluster
    requiredCols <- c(
        "avgLogFC",     # Seurat v2.1
        "cluster",      # Unmodified
        "description",  # Ensembl annotations
        "geneBiotype",  # Ensembl annotations
        "geneName",     # Renamed from `gene`
        "geneID",       # Ensembl annotations
        "pvalue"        # Renamed from `p_val`
    )
    assert_is_subset(requiredCols, colnames(markers))

    # Return as a tibble
    # Grouped by cluster
    # Arranged by P value (per cluster)
    markers %>%
        remove_rownames() %>%
        as_tibble() %>%
        # Ensure all the annotations added are camelCase
        camel(strict = FALSE) %>%
        # `padj` should come at the end, but isn't in legacy output
        # TODO Switch to base R method here
        select(
            c(
                "cluster",
                "geneID",
                "geneName",
                "pct1",
                "pct2",
                "avgLogFC",
                "pvalue"
            ),
            everything()
        ) %>%
        group_by(.data[["cluster"]]) %>%
        # Arrange by P value
        arrange(!!sym("pvalue"), .by_group = TRUE)
}



# Methods ======================================================================
#' @rdname sanitizeMarkers
#' @export
setMethod(
    "sanitizeMarkers",
    signature(
        object = "seurat",
        markers = "data.frame"
    ),
    .sanitizeMarkersSeurat
)
