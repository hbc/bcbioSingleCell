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
    sparse_counts <- readFileByExtension(matrix_file)
    if (pipeline == "bcbio") {
        colnames(sparse_counts) <-
            readFileByExtension(col_file) %>%
            snake
        rownames(sparse_counts) <-
            readFileByExtension(row_file)
    } else if (pipeline == "cellranger") {
        # Named `barcodes.tsv` but not actually tab delimited
        colnames(sparse_counts) <-
            readFileByExtension(
                col_file,
                col_names = "cellular_barcode",
                col_types = "c") %>%
            pull("cellular_barcode") %>%
            snake
        rownames(sparse_counts) <-
            readFileByExtension(
                row_file,
                col_names = c("ensgene", "symbol"),
                col_types = "cc") %>%
            pull("ensgene")
    }


    # Cellular barcode sanitization =====
    # Ensure cellular barcodes are snake_case
    colnames(sparse_counts) <- snake(colnames(sparse_counts))

    # CellRanger outputs unnecessary trailing `-1`.
    if (pipeline == "cellranger" &
        all(str_detect(colnames(sparse_counts), "_1$"))) {
        colnames(sparse_counts) <-
            str_replace(colnames(sparse_counts), "_1$", "")
    }

    # Reformat to `[acgt]{8}_[acgt]{8}` instead of `[acgt]{16}`
    if (all(str_detect(colnames(sparse_counts), "^[acgt]{16}$"))) {
        colnames(sparse_counts) <- colnames(sparse_counts) %>%
            str_replace("^([acgt]{8})([acgt]{8})$", "\\1_\\2")
    }

    # Add sample name
    # 8 nucleotides: inDrop, Chromium
    # 6 nucleotides: SureCell
    if (all(str_detect(colnames(sparse_counts), "^[acgt]{6,8}_"))) {
        colnames(sparse_counts) <- colnames(sparse_counts) %>%
            paste(snake(sample_name), ., sep = "_")
    } else {
        stop("Failed to add sample name")
    }

    # Return as dgCMatrix, for improved memory overhead
    as(sparse_counts, "dgCMatrix")
}



#' @rdname sparse_counts
.sparse_counts_tx2gene <- function(sparse_counts, tx2gene) {
    sparse_counts <- .strip_transcript_versions(sparse_counts)
    if (!all(rownames(sparse_counts) %in% rownames(tx2gene))) {
        warning("tx2gene missing transcripts present in sparse counts")
    }
    t2g <- tx2gene[rownames(sparse_counts), ] %>%
        set_rownames(rownames(sparse_counts))
    if (any(is.na(t2g[["enstxp"]]))) {
        t2g <- rownames_to_column(t2g)
        mapped <- filter(t2g, !is.na(.data[["ensgene"]]))
        unmapped <- filter(t2g, is.na(.data[["ensgene"]])) %>%
            mutate(enstxp = .data[["rowname"]],
                   ensgene = .data[["rowname"]])
        warning(paste("Unmapped transcripts:",
                      toString(unmapped[["ensgene"]])))
        t2g <- bind_rows(mapped, unmapped) %>%
            column_to_rownames %>%
            .[rownames(sparse_counts), ]
    }
    if (!identical(rownames(sparse_counts), rownames(t2g))) {
        stop("tx2gene rowname mismatch")
    }
    message("Converting transcript-level counts to gene-level")
    sparse_counts %>%
        set_rownames(t2g[["ensgene"]]) %>%
        aggregate.Matrix(rownames(.), fun = "sum")
}



#' @rdname sparse_counts
#' @description Strip transcript versions.
#' @return Sparse counts matrix with the transcript versions stripped.
.strip_transcript_versions <- function(sparse_counts) {
    transcripts <- rownames(sparse_counts)
    # Pattern matching against Ensembl transcript IDs
    # http://www.ensembl.org/info/genome/stable_ids/index.html
    # Examples: ENST (human); ENSMUST (mouse)
    enstxp_pattern <- "^(ENS.*T\\d{11})\\.\\d+$"
    if (any(str_detect(transcripts, enstxp_pattern))) {
        transcripts <- str_replace(transcripts, enstxp_pattern, "\\1")
        rownames(sparse_counts) <- transcripts
    }
    if (any(str_detect(transcripts, "\\.\\d+$"))) {
        stop("Incomplete transcript version removal")
    }
    sparse_counts
}
