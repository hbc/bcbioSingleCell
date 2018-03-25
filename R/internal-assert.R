#' @importFrom assertive assert_all_are_dirs
#' @importFrom assertive assert_all_are_existing_files
#' @importFrom assertive assert_all_are_greater_than_or_equal_to
#' @importFrom assertive assert_all_are_in_left_open_range
#' @importFrom assertive assert_all_are_matching_regex
#' @importFrom assertive assert_all_are_non_negative
#' @importFrom assertive assert_all_are_positive
#' @importFrom assertive assert_any_are_matching_regex
#' @importFrom assertive assert_are_disjoint_sets
#' @importFrom assertive assert_are_intersecting_sets
#' @importFrom assertive assert_are_identical
#' @importFrom assertive assert_has_dimnames
#' @importFrom assertive assert_has_no_duplicates
#' @importFrom assertive assert_has_names
#' @importFrom assertive assert_has_rows
#' @importFrom assertive assert_is_a_bool
#' @importFrom assertive assert_is_a_number
#' @importFrom assertive assert_is_a_string
#' @importFrom assertive assert_is_all_of
#' @importFrom assertive assert_is_an_integer
#' @importFrom assertive assert_is_any_of
#' @importFrom assertive assert_is_character
#' @importFrom assertive assert_is_data.frame
#' @importFrom assertive assert_is_environment
#' @importFrom assertive assert_is_factor
#' @importFrom assertive assert_is_function
#' @importFrom assertive assert_is_list
#' @importFrom assertive assert_is_non_empty
#' @importFrom assertive assert_is_of_length
#' @importFrom assertive assert_is_subset
#' @importFrom assertive assert_is_tbl_df
#' @importFrom assertive has_dims
#' @importFrom assertive has_names
#' @importFrom assertive is_a_string
#' @importFrom assertive is_character
#'
#' @importFrom basejump assertHasRownames
#' @importFrom basejump assertIsAHeaderLevel
#' @importFrom basejump assertIsAStringOrNULL
#' @importFrom basejump assertIsAnImplicitInteger
#' @importFrom basejump assertIsAnImplicitIntegerOrNULL
#' @importFrom basejump assertIsCharacterOrNULL
#' @importFrom basejump assertIsColorScaleContinuousOrNULL
#' @importFrom basejump assertIsColorScaleDiscreteOrNULL
#' @importFrom basejump assertIsFillScaleDiscreteOrNULL
#' @importFrom basejump assertIsGene2symbol
#' @importFrom basejump assertIsTx2gene
#'
#' @importFrom bcbioBase assertFormalInterestingGroups
NULL



.isAggregate <- function(object, stop = FALSE) {
    logical <- "sampleNameAggregate" %in% colnames(object)
    if (
        identical(logical, FALSE) &&
        identical(stop, TRUE)
    ) {
        abort("`sampleNameAggregate` column is required")
    }
    logical
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
.isSanitizedMarkers <- function(
    object,
    package = "Seurat"
) {
    package <- match.arg(package)

    # General checks ===========================================================
    if (!is(object, "grouped_df")) {
        return(FALSE)
    }
    else if (
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
            "gene",       # Gene symbol, we'll rename to "geneName"
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
