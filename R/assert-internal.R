.assertHasIdent <- function(object) {
    assert_is_subset("ident", colnames(colData(object)))
}



.isAggregate <- function(object, stop = FALSE) {
    logical <- "aggregate" %in% colnames(object)
    if (
        identical(logical, FALSE) &&
        identical(stop, TRUE)
    ) {
        stop("`aggregate` column is required")
    }
    logical
}



.isSanitizedMarkers <- function(
    object,
    package = "Seurat"
) {
    package <- match.arg(package)

    # General checks ===========================================================
    if (!is(object, "grouped_df")) {
        return(FALSE)
    } else if (
        is.null(attr(object, "vars")) ||
        attr(object, "vars") != "cluster"
    ) {
        return(FALSE)
    } else if (!"geneID" %in% colnames(object)) {
        return(FALSE)
    }

    # Package-specific checks ==================================================
    if (package == "Seurat") {
        # Check for `Seurat::FindAllMarkers()` return.
        # These columns are output in an inconsistent format, so we'll sanitize
        # into lowerCamelCase.
        seuratBlacklist <- c(
            "avg_diff",   # Legacy, now "avg_logFC"
            "avg_logFC",  # Renamed in v2.1
            "gene",
            "p_val",      # We'll rename to pvalue, matching DESeq2
            "p_val_adj",  # New in v2.1, we'll rename to padj, matching DESeq2
            "pct.1",
            "pct.2"
        )
        if (any(seuratBlacklist %in% colnames(object))) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    }
}



# Figure out which gene ID type is in use (geneID or geneName)
# Return TRUE if we're working with symbols (not recommended but used by Seurat)
.useGeneName <- function(object) {
    geneName <- as.character(gene2symbol(object)[["geneName"]])
    assert_is_non_empty(geneName)
    any(geneName %in% rownames(object))
}
