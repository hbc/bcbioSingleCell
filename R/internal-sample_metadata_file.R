#' Prepare sample metadata from file
#'
#' Appends reverse complement sequences by matching against the `sequence`
#' column.
#'
#' @rdname sample_metadata_file
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param file Sample metadata file.
#' @param sample_dirs Sample directory paths.
#' @param pipeline Pipeline used to process the samples.
#' @param unique_names Ensure sample names are unique.
#'
#' @return [tibble].
.sample_metadata_file <- function(
    file,
    sample_dirs,
    pipeline,
    unique_names = FALSE) {
    meta <- read_file_by_extension(file) %>% snake

    # Check file integrity
    if (pipeline == "bcbio") {
        required_cols <- c("index", "sequence", "sample_name")
    } else if (pipeline == "cellranger") {
        required_cols <- "sample_id"
    } else {
        stop("Unsupported pipeline")
    }
    if (!all(required_cols %in% colnames(meta))) {
        stop(paste("Metadata file missing columns:",
                   toString(required_cols)))
    }

    # bcbio-nextgen sample_id construction
    if (pipeline == "bcbio") {
        meta <- meta %>%
            mutate(revcomp = vapply(.data[["sequence"]],
                                    revcomp,
                                    character(1L)))
        # Map sample names to revcomp barcode indexes
        if ("file_name" %in% colnames(meta)) {
            meta <- meta %>%
                mutate(sample_id = paste(.data[["file_name"]],
                                         .data[["revcomp"]],
                                         sep = "_"))
        } else {
            revcomp_matches <- pull(meta, "revcomp")
            sample_matches <- basename(sample_dirs) %>%
                str_match("(.*)-([ACGT]+)$") %>%
                as("tibble") %>%
                set_colnames(c("sample_id", "file_name", "revcomp")) %>%
                filter(.data[["revcomp"]] %in% syms(revcomp_matches))
            meta <- suppressWarnings(
                full_join(meta, sample_matches, by = "revcomp"))
        }
    }

    # Ensure `sample_id` is valid name then arrange
    meta <- meta %>%
        # Sanitize sample ID into snake_case
        mutate(sample_id = snake(.data[["sample_id"]])) %>%
        arrange(!!sym("sample_id"))

    # Check that `sample_id` matches `sample_dirs`
    if (!all(meta[["sample_id"]] %in% names(sample_dirs))) {
        stop("Sample directory names don't match the sample metadata file")
    }

    # Set required columns, if empty
    if (is.null(meta[["file_name"]])) {
        meta[["file_name"]] <- "all_samples"
    }
    if (is.null(meta[["sample_name"]])) {
        meta[["sample_name"]] <- meta[["sample_id"]]
    }

    if (isTRUE(unique_names)) {
        # Ensure unique sample names
        if (any(duplicated(meta[["sample_name"]]))) {
            meta[["sample_name"]] <-
                str_c(meta[["sample_name"]], " (", meta[["file_name"]], ")")
            if (any(duplicated(meta[["sample_name"]]))) {
                stop("Unique sample name generation failed")
            }
        }
    }

    # Return
    if (pipeline == "bcbio") {
        meta %>%
            tidy_select(c(meta_priority_cols,
                          "index",
                          "sequence",
                          "revcomp"),
                        everything())
    } else {
        meta %>%
            tidy_select(meta_priority_cols, everything())
    }
}
