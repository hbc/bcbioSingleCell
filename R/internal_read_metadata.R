#' Read metadata
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param file Metadata file. CSV and XLSX formats are supported.
#' @param pattern Apply grep pattern matching to samples
#' @param pattern_col Column in data frame used for pattern subsetting
#' @param lanes Number of lanes used to split the samples into technical
#'   replicates (\code{_LXXX}) suffix.
#'
#' @return Metadata data frame
#' @export
read_metadata <- function(
    file,
    pattern = NULL,
    pattern_col = "sample_name",
    lanes = NULL) {
    if (!file.exists(file)) {
        stop("File not found")
    }

    # Check for format, based on metadata input name
    if (str_detect(file, "indrop")) {
        type <- "indrop"
    } else if (str_detect(file, "dropseq")) {
        type <- "dropseq"
    } else if (str_detect(file, "seqwell")) {
        type <- "seqwell"
    } else if (str_detect(file, "chromium|tenx|10x")) {
        type <- "10x"
    } else {
        stop("Unknown platform, please rename metadata file")
    }
    message(paste(type, "format detected"))

    # Load XLSX or CSV dynamically
    if (grepl("\\.xlsx$", file)) {
        metadata <- read_excel(file)
    } else {
        metadata <- read_csv(file)
    }

    # First column must be the FASTQ file name
    names(metadata)[1] <- "file_name"

    metadata <- metadata %>%
        set_names_snake %>%
        .[!is.na(.$description), ] %>%
        .[order(.$description), ]

    # Join platform-specific metadata
    if (type == "indrop") {
        i5_counts <- indrop_i5_index_counts()
        if (!is.null(i5_counts)) {
            metadata <- left_join(metadata, i5_counts)
        }
    }

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

    # Convert to data frame and set rownames
    metadata <- metadata %>% as.data.frame %>% set_rownames(.$sample_barcode)

    return(metadata)
}
