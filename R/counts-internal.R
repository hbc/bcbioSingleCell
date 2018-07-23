# Read bcbio Sparse Counts
# Matrix Market Exchange (MEX/MTX) format
.readBcbioSparseCounts <- function(sampleDir) {
    assert_is_a_string(sampleDir)
    assert_all_are_dirs(sampleDir)

    sampleStem <- basename(sampleDir)
    matrixFile <- file.path(sampleDir, paste0(sampleStem, ".mtx"))
    rowFile <- paste0(matrixFile, ".rownames")  # features
    colFile <- paste0(matrixFile, ".colnames")  # barcodes
    assert_all_are_existing_files(c(matrixFile, rowFile, colFile))

    # Attempt to load the column and rowname files first. If they're empty,
    # skip loading the MatrixMarket file, which will error otherwise. The bcbio
    # pipeline will output empty files for very low quality samples with no
    # cells that pass filtering.
    rownames <- read_lines(rowFile)
    colnames <- read_lines(colFile)
    if (!length(rownames) || !length(colnames)) {
        return(NULL)
    }

    counts <- readMM(matrixFile)

    assert_are_identical(length(rownames), nrow(counts))
    assert_are_identical(length(colnames), ncol(counts))
    rownames(counts) <- rownames
    colnames(counts) <- colnames

    counts
}



# Read Cell Ranger HDF5 Counts
# @seealso `cellrangerRkit::get_matrix_from_h5()`
.readCellRangerHDF5Counts <- function(sampleDir, filtered = TRUE) {
    assert_all_are_dirs(sampleDir)
    assert_is_a_string(sampleDir)
    assert_is_a_bool(filtered)

    if (isTRUE(filtered)) {
        prefix <- "filtered"
    } else {
        prefix <- "raw"
    }

    file <- file.path(
        sampleDir,
        "outs",
        paste0(prefix, "_gene_bc_matrices_h5.h5")
    )
    assert_all_are_existing_files(file)
    message(file)

    genomeBuild <- names(h5dump(file, load = FALSE))
    assert_is_a_string(genomeBuild)

    # name e.g. "/GRCh38/data"
    h5 <- h5read(file = file, name = genomeBuild)

    # CsparseMatrix instead of TsparseMatrix
    counts <- sparseMatrix(
        i = h5[["indices"]] + 1L,
        p = h5[["indptr"]],
        x = as.numeric(h5[["data"]]),
        dims = h5[["shape"]],
        giveCsparse = TRUE
    )

    rownames <- h5[["genes"]]
    colnames <- h5[["barcodes"]]
    # Move the multiplexed sample index number to the beginning
    colnames <- sub("^([ACGT]+)-(\\d+)$", "\\2-\\1", colnames)
    assert_are_identical(length(rownames), nrow(counts))
    assert_are_identical(length(colnames), ncol(counts))
    rownames(counts) <- rownames
    colnames(counts) <- colnames

    counts
}



# Read Cell Ranger Sparse Counts
# Matrix Market Exchange (MEX/MTX) format
.readCellRangerSparseCounts <- function(sampleDir, filtered = TRUE) {
    assert_is_a_string(sampleDir)
    assert_all_are_dirs(sampleDir)
    assert_is_a_bool(filtered)

    if (isTRUE(filtered)) {
        prefix <- "filtered"
    } else {
        prefix <- "raw"
    }

    countsDir <- file.path(
        sampleDir,
        "outs",
        paste0(prefix, "_gene_bc_matrices")
    )
    message(countsDir)
    assert_all_are_dirs(countsDir)

    genomeBuild <- list.dirs(
        path = countsDir,
        full.names = FALSE,
        recursive = FALSE
    )
    assert_is_a_string(genomeBuild)

    matrixFile <- file.path(countsDir, genomeBuild, "matrix.mtx")
    rowFile <- file.path(dirname(matrixFile), "genes.tsv")
    colFile <- file.path(dirname(matrixFile), "barcodes.tsv")
    assert_all_are_existing_files(c(matrixFile, rowFile, colFile))

    counts <- readMM(matrixFile)

    # `genes.tsv` is tab delimited
    rownames <- read_tsv(
        file = rowFile,
        col_names = c("geneID", "geneName"),
        col_types = "cc"
    )
    assert_has_rows(rownames)
    rownames <- pull(rownames, "geneID")

    # `barcodes.tsv` is not tab delimited
    colnames <- read_lines(colFile)
    # Move the multiplexed sample index number to the beginning for easier
    # downstream grep matching
    colnames <- sub("^([ACGT]+)-(\\d+)$", "\\2-\\1", colnames)

    assert_are_identical(length(rownames), nrow(counts))
    assert_are_identical(length(colnames), ncol(counts))
    rownames(counts) <- rownames
    colnames(counts) <- colnames

    counts
}



# Main internal read counts function, supporting bcbio and Cell Ranger
.readCounts <- function(
    sampleDirs,
    pipeline = c("bcbio", "cellranger"),
    # bcbio currently only supports MatrixMarket
    format = c("mtx", "hdf5"),
    # This argument only applies to Cell Ranger
    filtered = TRUE
) {
    assert_all_are_dirs(sampleDirs)
    assert_has_names(sampleDirs)
    pipeline <- match.arg(pipeline)
    format <- match.arg(format)
    assert_is_a_bool(filtered)
    message(paste("Reading counts in", format, "format"))

    list <- mapply(
        sampleDir = sampleDirs,
        sampleID = names(sampleDirs),
        MoreArgs = list(
            pipeline = pipeline,
            format = format
        ),
        FUN = function(
            sampleDir,
            sampleID,
            pipeline,
            format
        ) {
            if (pipeline == "bcbio") {
                counts <- .readBcbioSparseCounts(sampleDir)
            } else if (pipeline == "cellranger" && format == "mtx") {
                counts <- .readCellRangerSparseCounts(
                    sampleDir = sampleDir,
                    filtered = filtered
                )
            } else if (pipeline == "cellranger" && format == "hdf5") {
                counts <- .readCellRangerHDF5Counts(
                    sampleDir = sampleDir,
                    filtered = filtered
                )
            }

            # Consistently output as CsparseMatrix
            counts <- as(counts, "dgCMatrix")

            # Prefix cell barcodes with sample ID
            colnames <- colnames(counts) %>%
                paste(sampleID, ., sep = "_") %>%
                makeNames(unique = TRUE)
            colnames(counts) <- colnames

            # Ensure rownames are valid
            rownames <- rownames(counts) %>%
                makeNames(unique = TRUE)
            rownames(counts) <- rownames

            counts
        },
        SIMPLIFY = FALSE,
        USE.NAMES = TRUE
    )
    # Remove any empty items in list, which can result from low quality samples
    # with empty matrices in bcbio pipeline
    list <- Filter(Negate(is.null), list)

    do.call(cbind, list)
}
