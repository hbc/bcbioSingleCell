#' Read metadata.
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#' @param pattern Apply grep pattern matching to samples.
#' @param pattern_col Column in data frame used for pattern subsetting.
#'
#' @return Metadata data frame.
read_metadata <- function(
    run,
    pattern = NULL,
    pattern_col = "sample_name") {
    # Don't use `check_run(run)` here because the metadata will be missing
    file <- run$metadata_file
    lanes <- run$lanes

    # Load XLSX or CSV dynamically
    if (grepl("\\.xlsx$", file)) {
        message("Reading metadata XLSX file...")
        metadata <- read_excel(file)
    } else if (grepl("\\.csv$", file)) {
        message("Reading metadata CSV file...")
        metadata <- read_csv(file)
    } else {
        stop("Unsupported metadata file type")
    }

    # Check for platform-specific metadata, based on file name. We can
    # improve this step in the future using a YAML-based method instead.
    if (str_detect(file, "indrop")) {
        type <- "inDrop"
    } else if (str_detect(file, "dropseq")) {
        type <- "Drop-seq"
    } else if (str_detect(file, "seqwell")) {
        type <- "Seq-Well"
    } else if (str_detect(file, "chromium|tenx|10x")) {
        type <- "10x Chromium"
    } else {
        stop("Unknown platform, please rename metadata file")
    }
    message(paste(type, "metadata detected"))

    # First column must be the FASTQ file name
    names(metadata)[1] <- "file_name"

    metadata <- metadata %>%
        set_names_snake %>%
        # Keep rows with a description
        .[!is.na(.$description), ]
        # Order by description
        # .[order(.$description), ]

    # Lane split, if desired
    if (is.numeric(lanes)) {
        lane <- paste0("L", str_pad(1:lanes, 3, pad = "0"))
        metadata <- metadata %>%
            group_by(!!sym("file_name")) %>%
            expand_(~lane) %>%
            left_join(metadata, by = "file_name") %>%
            ungroup %>%
            mutate(file_name = paste(.data$file_name, .data$lane, sep = "_"),
                   description = .data$file_name)
    }

    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        metadata <- metadata[str_detect(metadata[[pattern_col]], pattern), ]
    }

    # Reverse complement matches
    metadata$reverse_complement <- sapply(metadata$sequence, revcomp)

    # Sample barcode identifier
    metadata$sample_barcode <- paste(
        metadata$description, metadata$reverse_complement, sep = "-")

    # Join platform-specific metadata
    if (type == "inDrop") {
        # Join i5 index counts only on non-lanesplit samples
        if (is.null(lanes)) {
            i5_counts <- indrop_i5_index_counts()
            if (!is.null(i5_counts)) {
                message("inDrop i5 index counts log detected")
                metadata <- left_join(metadata,
                                      i5_counts,
                                      by = "reverse_complement")
            }
        }
    }

    # Convert to data frame and set rownames
    metadata <- metadata %>%
        as.data.frame %>%
        set_rownames(.$sample_barcode)

    return(metadata)
}
