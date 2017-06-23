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
#' @note bcbio-nextgen outputs read counts at transcript level.
#'
#' @rdname sparse_counts
#' @keywords internal
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param matrix_file MatrixMart file to read.
#'
#' @return Sparse counts matrix.
#'
#' @seealso `help("dgCMatrix-class")`
.sparse_counts <- function(matrix_file) {
    name <- deparse(substitute(matrix_file))

    if (!file.exists(matrix_file)) {
        stop("Count matrix file missing")
    }

    # Detect a compressed MatrixMarket file
    if (str_detect(matrix_file, "\\.(bz2|gz|xz)$")) {
        # TODO Automatic compressed file handling
        stop("Compressed matrix handling will be added in a future update")
        # New draft method: pipe a connection to readMM
        # gzfile, bzfile, xzfile
        # gzcon, open(con, "r")

        # Old method: decompress on disk
        # matrix_file <- gunzip(matrix_file)
    }

    # Check for pipeline by support file structure
    parent_dir <- dirname(matrix_file)
    if (all(file.exists(paste0(matrix_file, ".colnames")) &
            file.exists(paste0(matrix_file, ".rownames")))) {
        pipeline <- "bcbio-nextgen"
        # Transcript-level counts
        count_level <- "transcript"
        col_file <- paste0(matrix_file, ".colnames")  # barcodes
        row_file <- paste0(matrix_file, ".rownames")  # transcripts
    } else if (all(file.exists(file.path(parent_dir, "barcodes.tsv")) &
                   file.exists(file.path(parent_dir, "genes.tsv")))) {
        # 10X Chromium CellRanger
        pipeline <- "cellranger"
        # Gene-level counts
        count_level <- "gene"
        col_file <- "barcodes.tsv"
        row_file <- "genes.tsv"
    } else {
        stop("Automatic pipeline detection failed")
    }

    message(paste(
        "Reading",
        name,
        paste0(count_level, "-level"),
        "counts from",
        pipeline))

    # Read the MatrixMarket file. Column names are molecular identifiers. Row
    # names are gene/transcript identifiers (depending on platform).
    counts <- readMM(matrix_file)
    rownames(counts) <- read_lines(row_file)
    colnames(counts) <- read_lines(col_file)

    # Count matrix sanitization
    if (pipeline == "bcbio-nextgen") {
        # inDrop v3 barcode sanitization. If cellular barcode filtering is
        # `auto`, the identifiers are output as `[ACGT]{16}` instead of
        # `[ACGT]{8}-[ACGT]{8}`. Check for this and correct, for consistent
        # barcode names. This should be fixed on the bcbio-nextgen pipline side
        # in a future update.
        if (any(str_detect(colnames(counts), "\\:[ACGT]{16}$"))) {
            colnames(counts) <- colnames(counts) %>%
                str_replace("\\:([ACGT]{8})([ACGT]{8})$", "\\:\\1-\\2")
        }
    } else if (pipeline == "cellranger") {
        # CellRanger outputs unnecessary trailing `-1` to barcodes.
        pattern <- "\\-1$"
        if (all(grepl(pattern, colnames(counts)))) {
            colnames(counts) <- str_replace(colnames(counts), pattern, "")
        }
    }

    # Return as dgCMatrix for improved memory overhead
    as(counts, "dgCMatrix")
}
