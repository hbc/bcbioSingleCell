# Read sample cellular barcode summary file
#
# @author Rory Kirchner
#
# @param file_name File name.
.barcode_file <- function(file_name) {
    df <- read_tsv(file_name,
                   col_names = c("cellular_barcode", "reads"),
                   col_types = cols(),
                   progress = FALSE)
    df$reads %>%
        as.numeric %>%
        set_names(df$cellular_barcode) %>%
        sort(decreasing = TRUE)
}



# Read [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen) output
#
# @author Michael Steinbaugh
# @author Rory Kirchner
#
# @param run [bcbioSCDataSet].
.barcodes <- function(run) {
    files <- run$sample_dirs %>%
        file.path(., paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        set_names(names(run$sample_dirs))
    if (!all(file.exists(files))) {
        warning("Barcode TSV file missing")
        return(NULL)
    }
    message("Reading barcode distributions...")
    pbmclapply(seq_along(files), function(a) {
        read_barcode_file(files[a])
    }) %>% set_names(names(files))
}



# Import transcript-level count data from a bcbio run into a sparse matrix
#
# @author Michael Steinbaugh
#
# @param strip_version Strip transcript version from identifier.
#
# @return Sparse counts matrix.
.counts <- function(run, strip_version = TRUE) {
    counts <- file.path(run$project_dir, "tagcounts.mtx") %>% read_counts
    if (isTRUE(strip_version)) {
        counts <- strip_transcript_versions(counts)
    }
    # Convert transcript-level counts to gene-level
    ensembl <- run$ensembl
    # Avoid setting rownames on data frames (tidy principle)
    indexes <- match(rownames(counts), ensembl$ensembl_transcript_id)
    gene_names <- ensembl %>%
        as.data.frame %>%
        .[indexes, "external_gene_name"]
    aggregate_sparse_features(counts, gene_names)
}







.program_versions <- function(project_dir) {
    file.path(project_dir, "programs.txt") %>%
        read_delim(",",
                   col_names = c("program", "version"),
                   col_types = "cc")
}



# Read a CSV file with readr, setting a given column as the rownames
#
# @author Rory Kirchner
#
# @param filename CSV to read.
# @param column Column to make into the rownames.
#
# @return Data frame.
.csv_with_rownames <- function(filename, column) {
    dat <- read_csv(filename, col_types = cols(), progress = FALSE) %>%
        as.data.frame
    rownames(dat) <- dat[, column]
    dat[, column] <- NULL
    dat
}



# Read metadata
#
# @author Michael Steinbaugh
#
# @param run bcbio-nextgen run.
# @param pattern Apply grep pattern matching to samples.
# @param pattern_col Column in data frame used for pattern subsetting.
#
# @return Metadata data frame.
.metadata <- function(
    run,
    pattern = NULL,
    pattern_col = "sample_name") {
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
        snake %>%
        filter(!is.na(.data$description))

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
                metadata <- left_join(
                    metadata, i5_counts, by = "reverse_complement")
            }
        }
    }

    # Convert to data frame and set rownames
    metadata %>%
        as.data.frame %>%
        set_rownames(.$sample_barcode)
}
