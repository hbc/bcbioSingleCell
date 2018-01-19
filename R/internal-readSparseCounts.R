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
#' @importFrom dplyr pull
#' @importFrom Matrix readMM
#' @importFrom readr read_lines read_tsv
#'
#' @author Michael Steinbaugh
#'
#' @param sampleDir Named character vector of sample directory containing the
#'   MatrixMart file.
#' @param pipeline Pipeline used to generate the MatrixMarket file. Defaults to
#'   bcbio-nextgen (`bcbio`). Also supports 10X Chromium CellRanger
#'   (`cellranger`).
#' @param umiType UMI type.
#'
#' @seealso `help("dgCMatrix-class")`
#'
#' @return `dgCMatrix`.
#' @noRd
.readSparseCounts <- function(
    sampleDir,
    pipeline,
    umiType) {
    sampleName <- names(sampleDir)
    if (is.null(sampleName)) {
        stop("Sample directory must be passed in as a named character vector",
             call. = FALSE)
    }
    message(sampleName)
    if (pipeline == "bcbio") {
        sampleStem <- basename(sampleDir)
        matrixFile <- file.path(sampleDir, paste0(sampleStem, ".mtx"))
        colFile <- paste0(matrixFile, ".colnames")  # barcodes
        rowFile <- paste0(matrixFile, ".rownames")  # transcripts
    } else if (pipeline == "cellranger") {
        filteredDir <- file.path(
            sampleDir,
            "outs",
            "filtered_gene_bc_matrices")
        matrixFile <- list.files(
            filteredDir,
            pattern = "matrix.mtx",
            full.names = TRUE,
            recursive = TRUE)
        colFile <- dirname(matrixFile) %>%
            file.path("barcodes.tsv")
        rowFile <- dirname(matrixFile) %>%
            file.path("genes.tsv")
    } else {
        stop("Unsupported pipeline", call. = FALSE)
    }

    # Check that all files exist
    if (!all(file.exists(matrixFile, colFile, rowFile))) {
        stop("Missing MatrixMarket file", call. = FALSE)
    }

    # Read the MatrixMarket file. Column names are molecular identifiers. Row
    # names are gene/transcript identifiers (depending on pipeline).
    counts <- readMM(matrixFile) %>%
        # Ensure dgCMatrix, for improved memory overhead
        as("dgCMatrix")

    if (pipeline == "bcbio") {
        colnames <- read_lines(colFile)
        rownames <- read_lines(rowFile)
    } else if (pipeline == "cellranger") {
        # Named `barcodes.tsv` but not actually tab delimited
        colnames <- read_tsv(
                file = colFile,
                col_names = "barcode",
                col_types = "c") %>%
            pull("barcode")
        rownames <- read_tsv(
                file = rowFile,
                col_names = c("ensgene", "symbol"),
                col_types = "cc") %>%
            pull("ensgene")
    }

    # Integrity checks
    if (!identical(length(colnames), ncol(counts))) {
        stop("Barcodes file doesn't match counts matrix rows")
    }
    if (!identical(length(rownames), nrow(counts))) {
        stop("Genes file doesn't match counts matrix columns")
    }

    colnames(counts) <- colnames %>%
        # Append sample name
        paste(sampleName, ., sep = "_") %>%
        # Convert dashes to underscores
        gsub(x = .,
             pattern = "-",
             replacement = "_") %>%
        make.names(unique = TRUE)
    rownames(counts) <- rownames %>%
        # Convert dashes to underscores
        gsub(x = .,
             pattern = "-",
             replacement = "_") %>%
        make.names(unique = TRUE)

    counts
}
