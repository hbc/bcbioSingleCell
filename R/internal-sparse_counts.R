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
#' @param matrix_file MatrixMart file to read.
#' @param pipeline Pipeline used to generate the MatrixMarket file. Defaults to
#'   bcbio-nextgen (`bcbio`). Also supports 10X Chromium CellRanger
#'   (`cellranger`).
#'
#' @return Sparse counts matrix.
#'
#' @seealso `help("dgCMatrix-class")`
.sparse_counts <- function(matrix_file, pipeline = "bcbio") {
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

    parent_dir <- dirname(matrix_file)

    if (pipeline == "bcbio") {
        count_level <- "transcript"
        col_file <- paste0(matrix_file, ".colnames")  # barcodes
        row_file <- paste0(matrix_file, ".rownames")  # transcripts
    } else if (pipeline == "cellranger") {
        count_level <- "gene"
        col_file <- file.path(parent_dir, "barcodes.tsv")
        row_file <- file.path(parent_dir, "genes.tsv")
    } else {
        stop("Unknown pipeline")
    }

    message(paste("Reading", name,
                  paste0(count_level, "-level"),
                  "counts from", pipeline))

    # Read the MatrixMarket file. Column names are molecular identifiers. Row
    # names are gene/transcript identifiers (depending on pipeline).
    counts <- readMM(matrix_file)
    rownames(counts) <- read_lines(row_file)
    colnames(counts) <- read_lines(col_file)

    # Count matrix sanitization
    if (pipeline == "bcbio") {
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
