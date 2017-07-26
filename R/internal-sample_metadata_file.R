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
    pipeline) {
    meta <- readFileByExtension(file) %>% snake

    if (pipeline == "bcbio") {
        # Rename legacy `samplename` column, if set
        if ("samplename" %in% colnames(meta)) {
            meta <- rename(meta, file_name = .data[["samplename"]])
        }

        # Rename `description` to `sample_name`, if set.
        if ("description" %in% colnames(meta)) {
            meta <- rename(meta, sample_name = .data[["description"]])
        }

        # Check for general required columns
        if (!all(c("file_name", "sample_name") %in% colnames(meta))) {
            stop("`file_name` and `sample_name` are required")
        }

        # Check if samples are demultiplexed
        if (length(unique(meta[["file_name"]])) == nrow(meta)) {
            # SureCell
            demultiplexed <- TRUE
        } else {
            # inDrop
            demultiplexed <- FALSE
        }

        if (isTRUE(demultiplexed)) {
            meta[["sample_id"]] <- meta[["sample_name"]]
        } else {
            # Match the UMI demultiplexed sample directories (e.g. inDrop)
            if (!"sequence" %in% colnames(meta)) {
                stop("Index i7 sequence required to generate `sample_id`")
            }
            meta <- meta %>%
                mutate(revcomp = vapply(.data[["sequence"]],
                                        revcomp,
                                        character(1L)),
                       sample_id = paste(.data[["file_name"]],
                                         .data[["revcomp"]],
                                         sep = "_"))
        }
    } else if (pipeline == "cellranger") {
        if (!"file_name" %in% colnames(meta)) {
            stop("Required `file_name` column missing")
        }
        meta[["sample_id"]] <- meta[["file_name"]]
    } else {
        stop("Unsupported pipeline")
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

    # Return
    tidy_select(meta, meta_priority_cols, everything())
}
