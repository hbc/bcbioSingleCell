#' Check for Sanitized Markers Data
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param object Input should be a sanitized markers [data.frame], returned from
#'   [sanitizeMarkers()].
#' @param package Name of package used to generate the markers.
#'
#' @return `TRUE` if sanitized, `FALSE` if not.
.checkSanitizedMarkers <- function(
    object,
    package = "Seurat") {
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
            FALSE
        } else {
            # Sanitized
            TRUE
        }
    } else {
        stop("Unsupported package")
    }
}
