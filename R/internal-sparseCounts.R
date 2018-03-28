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

    # Attempt to load the column and rowname files first. If they're empty,
    # skip loading the MatrixMarket file, which will error otherwise. The bcbio
    # pipeline outputs empty files.
    if (pipeline == "bcbio") {
        colnames <- read_lines(colFile)
        rownames <- read_lines(rowFile)
    } else if (pipeline == "cellranger") {
        # `barcodes.tsv` is not tab delimited
        colnames <- read_lines(colFile)

        # `genes.tsv` is tab delimited
        rownames <- read_tsv(
                file = rowFile,
                col_names = c("geneID", "geneName"),
                col_types = "cc"
            )
        assert_has_rows(rownames)
        rownames <- pull(rownames, "geneID")
    }

    # Early return on empty colnames
    if (!length(colnames)) {
        warn(paste(
            sampleID, "does not contain any cells"
        ))
        return(NULL)
    }

    # Read the MatrixMarket file.
    # Column names are molecular identifiers.
    # Row names are gene/transcript identifiers.
    # Always return dgCMatrix, for consistency and improved memory overhead.
    counts <- readMM(matrixFile) %>%
        as("dgCMatrix")

    # Integrity checks
    assert_are_identical(length(colnames), ncol(counts))
    assert_are_identical(length(rownames), nrow(counts))

    # Append `sampleID` to colnames and make valid
    colnames(counts) <- colnames %>%
        paste(sampleID, ., sep = "_") %>%
        makeNames(unique = TRUE)

    # Ensure rownames are valid
    rownames(counts) <- makeNames(rownames, unique = TRUE)

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
