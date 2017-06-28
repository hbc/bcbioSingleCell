#' Sample metadata
#'
#' Appends reverse complement sequences by matching against the `sequence`
#' column.
#'
#' @rdname sample_metadata
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param metadata Metadata data frame.
#' @param sample_dirs (*Optional*). Named character vector of sample directory
#'   paths.
#'
#' @return [tibble].
.sample_metadata <- function(file, sample_dirs) {
    # Read the metadata file
    metadata <- .read_file(file)
    required_cols <- c("index", "sequence", "sample_name")
    if (!all(required_cols %in% colnames(metadata))) {
        stop(paste("Metadata file missing required columns:",
                   toString(required_cols)))
    }
    metadata <- mutate(
        metadata,
        revcomp = vapply(.data[["sequence"]], revcomp, character(1L)))

    # Map sample names to inDrop barcode indexes (revcomp)
    if ("file_name" %in% colnames(metadata)) {
        # Sample barcode is structured as `file_name-revcomp`.
        metadata <- metadata %>%
            mutate(sample_id = paste(
                .data[["file_name"]], .data[["revcomp"]], sep = "-"))
    } else {
        # If file names aren't specified in the metadata, we assume the inDrop
        # indexes used are unique. Therefore, we can perform a full join against
        # the revcomp sequences present in the sample directories.
        revcomp_matches <- pull(metadata, "revcomp")
        sample_matches <- names(sample_dirs) %>%
            str_match("(.*)-([ACGT]+)$") %>%
            as_tibble %>%
            set_colnames(c("sample_id", "file_name", "revcomp")) %>%
            filter(.data[["revcomp"]] %in% syms(revcomp_matches))
        metadata <- suppressWarnings(
            full_join(metadata, sample_matches, by = "revcomp"))
    }

    metadata %>%
        tidy_select(c("file_name",
                 "index",
                 "sequence",
                 "revcomp",
                 "sample_id",
                 "sample_name"),
               everything()) %>%
        arrange(!!sym("sample_id"))
}
