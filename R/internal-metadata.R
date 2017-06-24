# FIXME Rory applied some updates to the MASTER branch...merge here
#' Prepare inDrop metadata
#'
#' Appends reverse complement sequences by matching against the `sequence`
#' column.
#'
#' @rdname indrop_metadata
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param metadata Metadata data frame.
#' @param sample_dirs Named character vector of sample directory paths.
#'
#' @return [DataFrame].
.indrop_metadata <- function(file, sample_dirs) {
    # Read the metadata file
    meta <- .read_file(file) %>% as.data.frame
    required_cols <- c("sample_name", "index", "sequence")
    if (!all(required_cols %in% colnames(meta))) {
        stop(paste("Metadata file missing required columns:",
                   toString(required_cols)))
    }
    meta <- mutate(
        meta,
        revcomp = vapply(.data[["sequence"]], revcomp, "character"))

    # Map sample names to inDrop barcode indexes (revcomp)
    if ("file_name" %in% colnames(meta)) {
        # Sample barcode is structured as `file_name-revcomp`.
        meta <- meta %>%
            mutate(sample_barcode = paste(.data[["file_name"]],
                                          .data[["revcomp"]],
                                          sep = "-"))
    } else {
        # If file names aren't specified in the metadata, we assume the inDrop
        # indexes used are unique. Therefore, we can perform a full join against
        # the revcomp sequences present in the sample directories.
        samples <- names(sample_dirs) %>%
            str_match("(.*)-([ACGT]{8})$") %>%
            as.data.frame %>%
            set_colnames(c("sample_barcode", "file_name", "revcomp"))
        meta <- full_join(meta, samples, by = "revcomp")
    }

    meta %>%
        tidy_select(c("file_name",
                 "index",
                 "sequence",
                 "revcomp",
                 "sample_barcode",
                 "sample_name"),
               everything()) %>%
        arrange(!!sym("sample_barcode")) %>%
        set_rownames(.[["sample_barcode"]]) %>%
        DataFrame
}
