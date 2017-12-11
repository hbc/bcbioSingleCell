#' Check for Sanitized Markers Data
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param object Input should be a sanitized markers [data.frame], returned from
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
            stop("Object must be 'grouped_df' class", call. = FALSE)
        } else {
            return(FALSE)
        }
    }
    if (is.null(attr(object, "vars")) |
        attr(object, "vars") != "cluster") {
        if (isTRUE(stop)) {
            stop("Object must be grouped by 'cluster'", call. = FALSE)
        } else {
            return(FALSE)
        }
    }
    if (!"ensgene" %in% colnames(object)) {
        if (isTRUE(stop)) {
            stop("Object missing 'ensgene' column", call. = FALSE)
        } else {
            return(FALSE)
        }
    }

    # Package-specific checks ==================================================
    if (package == "Seurat") {
        # Check for `Seurat::FindAllMarkers()` return.
        # These columns are output in an inconsistent format, so we'll sanitize
        # into lowerCamelCase.
        seuratBlacklist <-
            c("avg_diff",   # legacy, now "avg_logFC"
              "avg_logFC",  # renamed in v2.1
              "gene",       # gene symbol, we'll rename to "symbol"
              "p_val",      # we'll rename to pvalue, matching DESeq2
              "p_val_adj",  # new in v2.1, we'll rename to padj, matching DESeq2
              "pct.1",
              "pct.2")
        if (any(seuratBlacklist %in% colnames(object))) {
            # Original
            sanitized <- FALSE
        } else {
            # Sanitized
            sanitized <- TRUE
        }
    } else {
        stop("Unsupported package")
    }

    # Return or stop
    # Silent return if sanitized = TRUE and stop = FALSE
    if (identical(stop, FALSE)) {
        sanitized
    } else if (identical(stop, TRUE) &
               identical(sanitized, FALSE)) {
        stop(paste(
            "Markers table doesn't pass sanitization checks.",
            "Ensure that 'sanitizeMarkers()' has been run."
        ), call. = FALSE)
    }
}
