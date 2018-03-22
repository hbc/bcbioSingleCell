#' Read Sparse Counts
#'
#' Read single-cell RNA-seq counts saved as a
#' [MatrixMarket](https://people.sc.fsu.edu/~jburkardt/data/mm/mm.html) file
#' into a sparse counts matrix (`dgCMatrix`). Transcript ([rownames()]) and
#' molecular barcodes ([colnames()]) identifiers are set from the corresponding
#' dependency files automatically.
#'
#' This function supports automatic handling of compressed matrix files.
#'
#' @note bcbio-nextgen outputs counts at transcript level. 10X Chromium
#'   CellRanger outputs counts at gene level.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams general
#' @param sampleID Sample identifier. Must match the `sampleID` column of the
#'   sample metadata `data.frame`.
#' @param sampleDir Named character vector of sample directory containing the
#'   MatrixMart file.
#' @param umiType UMI type.
#'
#' @seealso `help("dgCMatrix-class")`
#'
#' @return `dgCMatrix`.
.readSparseCounts <- function(
    sampleID,
    sampleDir,
    pipeline = c("bcbio", "cellranger"),
    umiType
) {
    assert_is_a_string(sampleID)
    assert_is_a_string(sampleDir)
    assert_all_are_dirs(sampleDir)
    assert_is_a_string(pipeline)
    pipeline <- match.arg(pipeline)
    assert_is_a_string(umiType)

    if (pipeline == "bcbio") {
        sampleStem <- basename(sampleDir)
        matrixFile <- file.path(sampleDir, paste0(sampleStem, ".mtx"))
        colFile <- paste0(matrixFile, ".colnames")  # barcodes
        rowFile <- paste0(matrixFile, ".rownames")  # transcripts
    } else if (pipeline == "cellranger") {
        matrixFile <- list.files(
            sampleDir,
            pattern = "matrix.mtx",
            full.names = TRUE,
            recursive = TRUE
        )
        # Ensure we're using the filtered matrix
        matrixFile <- matrixFile[grepl(
            x = matrixFile,
            pattern = "filtered_gene_bc_matrices"
        )]
        if (length(matrixFile) != 1L) {
            abort("Failed to detect filtered matrix file")
        }
        colFile <- file.path(dirname(matrixFile), "barcodes.tsv")
        rowFile <- file.path(dirname(matrixFile), "genes.tsv")
    }

    # Check that all files exist
    assert_all_are_existing_files(c(matrixFile, colFile, rowFile))

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
                col_types = "c"
            ) %>%
            pull("barcode")
        rownames <- read_tsv(
                file = rowFile,
                col_names = c("geneID", "geneName"),
                col_types = "cc"
            ) %>%
            pull("geneID")
    }

    # Integrity checks
    assert_are_identical(length(colnames), ncol(counts))
    assert_are_identical(length(rownames), nrow(counts))

    colnames(counts) <- colnames %>%
        # Append sample name
        paste(sampleID, ., sep = "_") %>%
        # Convert dashes to underscores
        gsub(
            x = .,
            pattern = "-",
            replacement = "_"
        ) %>%
        make.names(unique = TRUE)
    rownames(counts) <- rownames %>%
        # Convert dashes to underscores
        gsub(
            x = .,
            pattern = "-",
            replacement = "_"
        ) %>%
        make.names(unique = TRUE)

    counts
}



.sparseCountsList <- function(sampleDirs, pipeline, umiType) {
    mcmapply(
        sampleID = names(sampleDirs),
        sampleDir = sampleDirs,
        MoreArgs = list(pipeline = pipeline, umiType = umiType),
        FUN = function(sampleID, sampleDir, pipeline, umiType) {
            .readSparseCounts(
                sampleID = sampleID,
                sampleDir = sampleDir,
                pipeline = pipeline,
                umiType = umiType
            )
        },
        SIMPLIFY = FALSE,
        USE.NAMES = TRUE
    )
}
