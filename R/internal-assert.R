# TODO Convert these functions into better assert checks



.checkAggregate <- function(object, stop = FALSE) {
    logical <- "sampleNameAggregate" %in% colnames(object)
    if (identical(logical, FALSE) &
        identical(stop, TRUE)) {
        abort("`sampleNameAggregate` column is required")
    }
    logical
}



.checkFilterParams <- function(object) {
    params <- metadata(object)[["filterParams"]]
    if (is.null(params)) {
        abort("`filterCells()` hasn't been applied to this dataset")
    }

    # `filterParams` metadata was stored as a named numeric vector up until
    # v0.0.28. We changed to storing as a list in v0.0.29, to allow for per
    # sample cutoffs.
    if (is.numeric(params)) {
        params <- as.list(params)
    }

    # Ensure all params are numeric
    if (!all(vapply(
        X = params,
        FUN = is.numeric,
        FUN.VALUE = logical(1L)
    ))) {
        abort("Filter parameters must be numeric")
    }
}



#' Check for Sanitized Markers Data
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param object Input should be a sanitized markers `data.frame`, returned from
#'   [sanitizeMarkers()].
#' @param package Name of package used to generate the markers.
#' @param stop Stop on failure, instead of returning `FALSE`.
#'
#' @return `TRUE` if sanitized, `FALSE` if not.
.checkSanitizedMarkers <- function(
    object,
    package = "Seurat",
    stop = FALSE) {
    # General checks ===========================================================
    if (!is(object, "grouped_df")) {
        if (isTRUE(stop)) {
            abort("Object must be `grouped_df` class")
        } else {
            return(FALSE)
        }
    }
    if (is.null(attr(object, "vars")) |
        attr(object, "vars") != "cluster") {
        if (isTRUE(stop)) {
            abort("Object must be grouped by `cluster` column")
        } else {
            return(FALSE)
        }
    }
    if (!"geneID" %in% colnames(object)) {
        if (isTRUE(stop)) {
            abort("Object missing `geneID` column")
        } else {
            return(FALSE)
        }
    }

    # Package-specific checks ==================================================
    if (package == "Seurat") {
        # Check for `Seurat::FindAllMarkers()` return.
        # These columns are output in an inconsistent format, so we'll sanitize
        # into lowerCamelCase.
        seuratBlacklist <- c(
            "avg_diff",   # Legacy, now "avg_logFC"
            "avg_logFC",  # Renamed in v2.1
            "gene",       # Gene symbol, we'll rename to "geneName"
            "p_val",      # We'll rename to pvalue, matching DESeq2
            "p_val_adj",  # New in v2.1, we'll rename to padj, matching DESeq2
            "pct.1",
            "pct.2"
        )
        if (any(seuratBlacklist %in% colnames(object))) {
            # Original
            sanitized <- FALSE
        } else {
            # Sanitized
            sanitized <- TRUE
        }
    } else {
        abort("Unsupported package")
    }

    # Return or stop
    # Silent return if sanitized = TRUE and stop = FALSE
    if (identical(stop, FALSE)) {
        sanitized
    } else if (identical(stop, TRUE) & identical(sanitized, FALSE)) {
        abort(paste(
            "Markers table doesn't pass sanitization checks.",
            "Ensure that `sanitizeMarkers()` has been run."
        ))
    }
}
