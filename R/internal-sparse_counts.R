#' Read MatrixMarket file
#'
#' Load single-cell RNA-seq counts saved as a [MatrixMarket
#' file](https://people.sc.fsu.edu/~jburkardt/data/mm/mm.html) into a sparse
#' counts matrix (`dgCMatrix`). Transcript ([rownames]) and molecular barcodes
#' ([colnames]) identifiers are set from the corresponding dependency files
#' automatically.
#'
#' This function supports automatic handling of compressed matrix files.
#'
#' @note bcbio-nextgen outputs counts at transcript level. 10X Chromium
#'   CellRanger outputs counts at gene level.
#'
#' @rdname sparse_counts
#' @keywords internal
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param sample_dir Named character vector of sample directory containing the
#'   MatrixMart file.
#' @param pipeline Pipeline used to generate the MatrixMarket file. Defaults to
#'   bcbio-nextgen (`bcbio`). Also supports 10X Chromium CellRanger
#'   (`cellranger`).
#'
#' @return Sparse counts matrix.
#'
#' @seealso `help("dgCMatrix-class")`
.sparse_counts <- function(
    sample_dir,
    pipeline = "bcbio") {
    sample_name <- names(sample_dir)
    if (is.null(sample_name)) {
        stop("Sample directory must be passed in as a named character vector")
    }
    message(sample_name)
    if (pipeline == "bcbio") {
        sample_stem <- basename(sample_dir)
        matrix_file <- file.path(sample_dir, str_c(sample_stem, ".mtx"))
        col_file <- str_c(matrix_file, ".colnames")  # barcodes
        row_file <- str_c(matrix_file, ".rownames")  # transcripts
    } else if (pipeline == "cellranger") {
        matrix_file <- file.path(sample_dir, "matrix.mtx")
        col_file <- file.path(sample_dir, "barcodes.tsv")
        row_file <- file.path(sample_dir, "genes.tsv")
    } else {
        stop("Unsupported pipeline")
    }

    # Read the MatrixMarket file. Column names are molecular identifiers. Row
    # names are gene/transcript identifiers (depending on pipeline).
    sparse_counts <- read_file_by_extension(matrix_file)
    if (pipeline == "bcbio") {
        colnames(sparse_counts) <-
            read_file_by_extension(col_file)
        rownames(sparse_counts) <-
            read_file_by_extension(row_file)
    } else if (pipeline == "cellranger") {
        # Named `barcodes.tsv` but not actually tab delimited
        colnames(sparse_counts) <-
            read_file_by_extension(
                col_file,
                col_names = "cellular_barcode",
                col_types = "c") %>%
            pull("cellular_barcode")
        rownames(sparse_counts) <-
            read_file_by_extension(
                row_file,
                col_names = c("ensgene", "symbol"),
                col_types = "cc") %>%
            pull("ensgene")
    }


    # Cellular barcode sanitization =====
    # CellRanger outputs unnecessary trailing `-1`.
    if (pipeline == "cellranger" &
        all(str_detect(colnames(sparse_counts), "-1$"))) {
        colnames(sparse_counts) <-
            str_replace(colnames(sparse_counts), "-1$", "")
    }

    # Reformat to `[ACGT]{8}-[ACGT]{8}` instead of `[ACGT]{16}`
    if (all(str_detect(colnames(sparse_counts), "^[ACGT]{16}$"))) {
        colnames(sparse_counts) <- colnames(sparse_counts) %>%
            str_replace("^([ACGT]{8})([ACGT]{8})$", "\\1_\\2")
    }

    # Add sample name
    if (all(str_detect(colnames(sparse_counts), "^[ACGT]{8}"))) {
        colnames(sparse_counts) <- colnames(sparse_counts) %>%
            paste(make.names(sample_name), ., sep = "_")
    }

    # Return as dgCMatrix, for improved memory overhead
    as(sparse_counts, "dgCMatrix")
}



#' @rdname sparse_counts
.sparse_counts_tx2gene <- function(tx_sparse_counts, genome_build) {
    tx_sparse_counts <- .strip_transcript_versions(tx_sparse_counts)
    tx2gene <- tx2gene(genome_build) %>%
        .[rownames(tx_sparse_counts), ]
    if (any(is.na(tx2gene[["enstxp"]]))) stop("Unmapped transcripts present")
    message("Converting transcript-level counts to gene-level")
    tx_sparse_counts %>%
        set_rownames(tx2gene[["ensgene"]]) %>%
        aggregate.Matrix(rownames(.), fun = "sum")
}



#' @rdname sparse_counts
#' @description Strip transcript versions.
#' @return Sparse counts matrix with the transcript versions stripped.
.strip_transcript_versions <- function(sparse_counts) {
    transcripts <- rownames(sparse_counts)
    # Tight pattern matching against Ensembl stable transcript IDs
    # http://www.ensembl.org/info/genome/stable_ids/index.html
    enstxp_pattern <- "^(ENS[A-Z]{3}T\\d{11})\\.\\d+$"
    if (any(str_detect(transcripts, enstxp_pattern))) {
        transcripts <- str_replace(transcripts, enstxp_pattern, "\\1")
        rownames(sparse_counts) <- transcripts
    }
    sparse_counts
}
