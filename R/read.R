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
##' read by extension
##'
##' Can read in XLSX/CSV/TSV files using the appropriate reader
##' @keywords internal
##' @param filename A filename to read in
##'
##' @return data frame of the file
##' @author Rory Kirchner
read_by_extension <- function(filename) {
    ext = tools::file_ext(filename)
    if (ext == "xlsx") {
        message("Reading XLSX file...")
        metadata <- read_excel(filename)
    } else if (ext == "csv") {
        message("Reading CSV file...")
        metadata <- read_csv(filename)
    } else if (ext == "tsv") {
        message("Reading TSV file...")
        metadata <- read_tsv(filename)
    }
    else {
        stop("Unsupported file type")
    }
}

#' Read metadata
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#'
#' @return Metadata data frame.
read_metadata <- function(run) {
    sample_metadata_file <- run$sample_metadata_file
    well_metadata_file <- run$well_metadata_file
    cells <- colnames(run$counts)
    samples <- stringr::str_split_fixed(cells, ":", 2)[,1]
    wells <- stringr::str_split_fixed(cells, ":", 2)[,2]
    metadata <- data.frame(sample_id=samples, well_id=wells, cell_id=cells)
    default_sampledata <- stringr::str_split_fixed(samples, "-", 2)

    if(!is.null(sample_metadata_file)) {
      metadata <- metadata %>%
        left_join(read_by_extension(sample_metadata_file) %>%
                mutate(sample_id=as.factor(sample_id)), by="sample_id")
    }
    if(!is.null(well_metadata_file)) {
      metadata <- metadata %>%
        left_join(read_by_extension(well_metadata_file) %>%
                mutate(sample_id=as.factor(sample_id)), by="well_id")
    }
    if(!"sample_name" %in% colnames(metadata)) {
      metadata$sample_name <- samples
    }
    if(!"sample_barcode" %in% colnames(metadata)) {
      metadata$sample_barcode <- default_sampledata[,2]
    }
    # Convert to data frame and set rownames
    metadata %>%
        as.data.frame %>%
        set_rownames(.$cell_id)
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
    barcodes_file <- file.path(data_dir, "barcodes.tsv")
    genes_file <- file.path(data_dir, "genes.tsv")
    matrix_file <- file.path(data_dir, "matrix.mtx")

    if (!file.exists(barcodes_file)) {
        stop("Barcodes file missing")
    }
    if (!file.exists(genes_file)) {
        stop("Gene names file missing")
    }
    if (!file.exists(matrix_file)) {
        stop("MatrixMart counts file missing")
    }

    # Read MatrixMart expression counts matrix
    counts <- readMM(matrix_file)

    # Assign barcodes to colnames
    barcodes <- readLines(barcodes_file)
    # Check for `-1` suffix in barcode names, indicative of 10X data
    pattern <- "\\-1$"
    if (all(grepl(pattern, barcodes))) {
        barcodes <- gsub(pattern, "", barcodes)
    }
    colnames(counts) <- barcodes

    # Assign gene names (symbols) to rownames
    genes <- read_tsv(genes_file,
                      col_names = c("ensembl_gene_id",
                                    "external_gene_name"))
    rownames(counts) <- make.unique(genes$external_gene_name)

    # Coerce dgTMatrix to dgCMatrix
    as(counts, "dgCMatrix")
}
