#' Read Sparse Counts
#'
#' Read single-cell RNA-seq counts saved as a
#' [MatrixMarket](https://people.sc.fsu.edu/~jburkardt/data/mm/mm.html) file
#' into a sparse counts matrix (`dgCMatrix`). Transcript ([rownames]) and
#' molecular barcodes ([colnames]) identifiers are set from the corresponding
#' dependency files automatically.
#'
#' This function supports automatic handling of compressed matrix files.
#'
#' @note bcbio-nextgen outputs counts at transcript level. 10X Chromium
#'   CellRanger outputs counts at gene level.
#'
#' @rdname internal-readSparseCounts
#' @keywords internal
#'
#' @param sampleDir Named character vector of sample directory containing the
#'   MatrixMart file.
#' @param pipeline Pipeline used to generate the MatrixMarket file. Defaults to
#'   bcbio-nextgen (`bcbio`). Also supports 10X Chromium CellRanger
#'   (`cellranger`).
#'
#' @return `dgCMatrix`.
#'
#' @seealso `help("dgCMatrix-class")`
.readSparseCounts <- function(
    sampleDir,
    pipeline = "bcbio") {
    sampleName <- names(sampleDir)
    if (is.null(sampleName)) {
        stop("Sample directory must be passed in as a named character vector")
    }
    message(sampleName)
    if (pipeline == "bcbio") {
        sampleStem <- basename(sampleDir)
        matrixFile <- file.path(sampleDir, paste0(sampleStem, ".mtx"))
        colFile <- paste0(matrixFile, ".colnames")  # barcodes
        rowFile <- paste0(matrixFile, ".rownames")  # transcripts
    } else if (pipeline == "cellranger") {
        filteredDir <- file.path(sampleDir,
                                 "outs",
                                 "filtered_gene_bc_matrices")
        matrixFile <- list.files(filteredDir,
                                 pattern = "matrix.mtx",
                                 full.names = TRUE,
                                 recursive = TRUE)
        colFile <- dirname(matrixFile) %>%
            file.path("barcodes.tsv")
        rowFile <- dirname(matrixFile) %>%
            file.path("genes.tsv")
    } else {
        stop("Unsupported pipeline")
    }

    # Check that all files exist
    if (!all(file.exists(matrixFile, colFile, rowFile))) {
        stop("Missing MatrixMarket file")
    }

    # Read the MatrixMarket file. Column names are molecular identifiers. Row
    # names are gene/transcript identifiers (depending on pipeline).
    sparseCounts <- readFileByExtension(matrixFile)
    if (pipeline == "bcbio") {
        colnames(sparseCounts) <-
            readFileByExtension(colFile) %>%
            str_replace_all("-", "_")
        rownames(sparseCounts) <-
            readFileByExtension(rowFile)
    } else if (pipeline == "cellranger") {
        # Named `barcodes.tsv` but not actually tab delimited
        colnames(sparseCounts) <-
            readFileByExtension(
                colFile,
                col_names = "cellularBarcode",
                col_types = "c") %>%
            pull("cellularBarcode") %>%
            str_replace_all("-", "_")
        rownames(sparseCounts) <-
            readFileByExtension(
                rowFile,
                col_names = c("ensgene", "symbol"),
                col_types = "cc") %>%
            pull("ensgene")
    }

    # Cellular barcode sanitization =====
    # CellRanger outputs unnecessary trailing `-1`.
    if (pipeline == "cellranger" &
        all(str_detect(colnames(sparseCounts), "_1$"))) {
        colnames(sparseCounts) <-
            str_replace(colnames(sparseCounts), "_1$", "")
    }

    # Reformat to `[ACGT]{8}_[ACGT]{8}` instead of `[ACGT]{16}`
    if (all(str_detect(colnames(sparseCounts), "^[ACGT]{16}$"))) {
        colnames(sparseCounts) <- colnames(sparseCounts) %>%
            str_replace("^([ACGT]{8})([ACGT]{8})$", "\\1_\\2")
    }

    # Add sample name
    # 8 nucleotides: inDrop, Chromium
    # 6 nucleotides: SureCell
    if (all(str_detect(colnames(sparseCounts), "^[ACGT]+_"))) {
        colnames(sparseCounts) <- colnames(sparseCounts) %>%
            paste(sampleName, ., sep = "_")
    } else {
        stop("Failed to add sample name")
    }

    # Return as dgCMatrix, for improved memory overhead
    as(sparseCounts, "dgCMatrix")
}



#' Transcript To Gene-Level Sparse Counts
#'
#' @rdname internal-sparseCountsTx2Gene
#' @keywords internal
#'
#' @param sparseCounts Sparse counts matrix (`dgCMatrix`).
#' @param tx2gene Transcript to gene identifier mappings.
#'
#' @return `dgCMatrix`.
.sparseCountsTx2Gene <- function(sparseCounts, tx2gene) {
    sparseCounts <- .stripTranscriptVersions(sparseCounts)
    if (!all(rownames(sparseCounts) %in% rownames(tx2gene))) {
        missing <- rownames(sparseCounts) %>%
            .[!. %in% rownames(tx2gene)]
        stop(paste(
            length(missing),
            "transcripts in matrix missing from tx2gene:",
            toString(head(missing)),
            "..."))
    }
    message("Converting transcript-level counts to gene-level")
    t2g <- tx2gene[rownames(sparseCounts), ]
    sparseCounts %>%
        set_rownames(t2g[["ensgene"]]) %>%
        aggregate.Matrix(rownames(.), fun = "sum")
}



#' Strip Transcript Versions
#'
#' @rdname internal-stripTranscriptVersions
#' @keywords internal
#'
#' @param sparseCounts Sparse counts matrix (`dgCMatrix`).
#'
#' @return `dgCMatrix`.
.stripTranscriptVersions <- function(sparseCounts) {
    transcripts <- rownames(sparseCounts)
    # Pattern matching against Ensembl transcript IDs
    # http://www.ensembl.org/info/genome/stable_ids/index.html
    # Examples: ENST (human); ENSMUST (mouse)
    enstxpPattern <- "^(ENS.*T\\d{11})\\.\\d+$"
    if (any(str_detect(transcripts, enstxpPattern))) {
        transcripts <- str_replace(transcripts, enstxpPattern, "\\1")
        rownames(sparseCounts) <- transcripts
    }
    if (any(str_detect(transcripts, "\\.\\d+$"))) {
        stop("Incomplete transcript version removal")
    }
    sparseCounts
}
