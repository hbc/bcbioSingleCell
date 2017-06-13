#' Read sample cellular barcode summary file
#'
#' @keywords internal
#' @author Rory Kirchner
#'
#' @param file_name File name.
#'
#' @export
read_barcode_file <- function(file_name) {
    df <- read_tsv(file_name,
                   col_names = c("cellular_barcode", "reads"),
                   progress = FALSE)
    df$reads %>%
        as.numeric %>%
        set_names(df$cellular_barcode) %>%
        sort(decreasing = TRUE)
}



#' Read [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen) output
#'
#' @rdname read_bcbio
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run [bcbioSCDataSet].
#'
#' @export
read_bcbio_barcodes <- function(run) {
    files <- run$sample_dirs %>%
        file.path(., paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        set_names(names(run$sample_dirs))
    if (!all(file.exists(files))) {
        warning("Barcode TSV file missing")
        return(NULL)
    }
    message("Reading barcode distributions...")
    list <- pbmclapply(seq_along(files), function(a) {
        read_barcode_file(files[a])
    }) %>% set_names(names(files))
}



#' @rdname read_bcbio
#' @description Import transcript-level count data from a bcbio run into a
#'   sparse matrix.
#' @param strip_version Strip transcript version from identifier.
#' @return Sparse counts matrix.
#' @export
read_bcbio_counts <- function(run, strip_version = TRUE) {
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



#' @rdname read_bcbio
#' @export
read_bcbio_programs <- function(run) {
    file.path(run$project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))
}



#' Read a MatrixMart file, setting the row and column names
#'
#' Note that for a bcbio run, this function will return transcript-level counts.
#'
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param matrix_file MatrixMart file to read.
#'
#' @return Sparse counts matrix.
#' @export
read_counts <- function(matrix_file) {
    # Detect gzip file
    pattern <- "\\.gz$"
    if (str_detect(matrix_file, pattern)) {
        matrix_file <- gunzip(matrix_file)
    }

    row_file <- paste0(matrix_file, ".rownames")
    col_file <- paste0(matrix_file, ".colnames")

    if (!file.exists(matrix_file)) {
        stop("MatrixMart counts file missing")
    }
    if (!file.exists(row_file)) {
        stop("Rownames file could not be found")
    }
    if (!file.exists(col_file)) {
        stop("Colnames file could not be found")
    }

    counts <- readMM(matrix_file)
    rownames(counts) <- read_lines(row_file)
    colnames(counts) <- read_lines(col_file)

    # [fix] Correct malformed celluar barcodes.
    # Need to update on the bcbio-nextgen platform side.
    if (any(str_detect(colnames(counts), "\\:[ACGT]{16}$"))) {
        colnames(counts) <- colnames(counts) %>%
            str_replace("\\:([ACGT]{8})([ACGT]{8})$", "\\:\\1-\\2")
    }

    # Coerce dgTMatrix to dgCMatrix
    as(counts, "dgCMatrix")
}



#' Read a CSV file with readr, setting a given column as the rownames
#'
#' @keywords internal
#' @author Rory Kirchner
#'
#' @param filename CSV to read.
#' @param column Column to make into the rownames.
#'
#' @return Data frame.
#' @export
read_csv_with_rownames <- function(filename, column) {
    dat <- read_csv(filename, progress = FALSE) %>%
        as.data.frame
    rownames(dat) <- dat[, column]
    dat[, column] <- NULL
    dat
}



#' Read metadata
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



#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) output
#'
#' Read into sparse matrix from barcodes.tsv, genes.tsv, and matrix.mtx files.
#'
#' @author Michael Steinbaugh
#'
#' @param data_dir Data directory
#'
#' @return Sparse counts matrix (dgCMatrix).
#' @export
read_10x <- function(data_dir) {
    # Recurse through specific data directory and identify sample subdirectories
    # by the presence of a `matrix.mtx` counts file.
    sample_dirs <- list.files(
        data_dir,
        full.names = TRUE,
        pattern = "matrix.mtx",
        recursive = TRUE) %>%
        normalizePath %>%
        dirname %>%
        set_names(basename(.))
    message(paste(length(sample_dirs), "samples detected"))
    message("Reading 10X Chromium samples...")
    lst <- pbmclapply(seq_along(sample_dirs), function(a) {
        name <- names(sample_dirs)[a]
        sample_dir <- sample_dirs[a]
        barcodes_file <- file.path(sample_dir, "barcodes.tsv")
        if (!file.exists(barcodes_file)) {
            stop("Barcodes file missing")
        }
        genes_file <- file.path(sample_dir, "genes.tsv")
        if (!file.exists(genes_file)) {
            stop("Gene names file missing")
        }
        matrix_file <- file.path(sample_dir, "matrix.mtx")
        if (!file.exists(matrix_file)) {
            stop("MatrixMart counts file missing")
        }

        # Read MatrixMart expression counts matrix
        counts <- readMM(matrix_file)

        # Assign barcodes to colnames
        barcodes <- read_lines(barcodes_file)
        # Check for `-1` suffix in barcode names, indicative of 10X data
        pattern <- "\\-1$"
        if (all(grepl(pattern, barcodes))) {
            barcodes <- str_replace(barcodes, pattern, "")
        }
        # Append the sample name to the barcodes. This will be used after
        # list is generated to [cBind()] all counts to a single sparse matrix.
        colnames(counts) <- paste(name, barcodes, sep = "-")

        # Assign gene names (symbols) to rownames
        genes <- read_tsv(genes_file,
                          col_names = c("ensembl_gene_id",
                                        "external_gene_name"))
        rownames(counts) <- make.unique(genes$external_gene_name)

        # Coerce dgTMatrix to dgCMatrix (sparse)
        as(counts, "dgCMatrix")
    }) %>% set_names(names(sample_dirs))
    # [fix] Any reason why we'd want to keep the files separate? Does this make
    # individual sample analysis easier?
    message("Aggregating samples to a single sparse matrix...")
    do.call(cBind, lst)
}
