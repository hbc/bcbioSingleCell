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
#' @param sample_metadata_file Sample metadata file.
#' @param sample_dirs Sample directory paths.
#' @param pipeline Pipeline used to process the samples.
#'
#' @return [tibble].
.sample_metadata_file <- function(
    sample_metadata_file,
    sample_dirs,
    pipeline) {
    sample_metadata <- read_file_by_extension(sample_metadata_file)

    # Check sample_metadata file integrity
    if (pipeline == "bcbio") {
        required_cols <- c("index", "sequence", "sample_name")
    } else if (pipeline == "cellranger") {
        required_cols <- "sample_id"
    } else {
        stop("Unsupported pipeline")
    }
    if (!all(required_cols %in% colnames(sample_metadata))) {
        stop(paste("Metadata file missing columns:",
                   toString(required_cols)))
    }

    # bcbio-nextgen sample_id construction
    if (pipeline == "bcbio") {
        sample_metadata <- mutate(
            sample_metadata,
            revcomp = vapply(.data[["sequence"]], revcomp, character(1L)))
        # Map sample names to revcomp barcode indexes
        if ("file_name" %in% colnames(sample_metadata)) {
            sample_metadata <- sample_metadata %>%
                mutate(sample_id = paste(
                    .data[["file_name"]], .data[["revcomp"]], sep = "-"))
        } else {
            revcomp_matches <- pull(sample_metadata, "revcomp")
            sample_matches <- basename(sample_dirs) %>%
                str_match("(.*)-([ACGT]+)$") %>%
                as_tibble %>%
                set_colnames(c("sample_id", "file_name", "revcomp")) %>%
                filter(.data[["revcomp"]] %in% syms(revcomp_matches))
            sample_metadata <- suppressWarnings(
                full_join(sample_metadata, sample_matches, by = "revcomp"))
        }
        sample_metadata <- sample_metadata %>%
            tidy_select(c("file_name",
                          "index",
                          "sequence",
                          "revcomp",
                          "sample_id",
                          "sample_name"),
                        everything())
    }

    # Check that sample_ids match samples
    if (!all(sample_metadata[["sample_id"]] %in% basename(sample_dirs))) {
        stop("Sample directory names don't match the sample metadata file")
    }

    sample_metadata %>%
        arrange(!!sym("sample_id"))
}
