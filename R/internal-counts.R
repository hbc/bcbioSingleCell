# Read bcbio Sparse Counts
# Matrix Market Exchange (MEX/MTX) format
.readBcbioSparseCounts <- function(sampleDir) {
    assert_is_a_string(sampleDir)
    assert_all_are_dirs(sampleDir)

    sampleStem <- basename(sampleDir)
    matrixFile <- file.path(sampleDir, paste0(sampleStem, ".mtx"))
    colFile <- paste0(matrixFile, ".colnames")  # barcodes
    rowFile <- paste0(matrixFile, ".rownames")  # transcripts
    assert_all_are_existing_files(c(matrixFile, colFile, rowFile))

    # Attempt to load the column and rowname files first. If they're empty,
    # skip loading the MatrixMarket file, which will error otherwise. The bcbio
    # pipeline will output empty files for very low quality samples with no
    # cells that pass filtering.
    colnames <- read_lines(colFile)
    rownames <- read_lines(rowFile)
    if (!length(colnames) || !length(rownames)) {
        warning(paste(sampleID, "does not contain any cells"))
        return(NULL)
    }

    counts <- readMM(matrixFile)
    assert_are_identical(length(colnames), ncol(counts))
    assert_are_identical(length(rownames), nrow(counts))
    colnames(counts) <- colnames
    rownames(counts) <- rownames
    counts
}



# Read Cell Ranger HDF5 Counts
# filtered: sample/outs/filtered_gene_bc_matrices_h5.h5
#      raw: sample/outs/raw_gene_bc_matrices_h5.h5
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

    rownames(counts) <- h5[["genes"]]
    colnames(counts) <- h5[["barcodes"]]
    counts
}



# Read Cell Ranger Sparse Counts
# Matrix Market Exchange (MEX/MTX) format
# filtered: sample/outs/filtered_gene_bc_matrices/GRCh38/matrix.mtx
#      raw: sample/outs/raw_gene_bc_matrices/GRCh38/matrix.mtx
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
    colFile <- file.path(dirname(matrixFile), "barcodes.tsv")
    rowFile <- file.path(dirname(matrixFile), "genes.tsv")
    assert_all_are_existing_files(c(matrixFile, colFile, rowFile))

    # `barcodes.tsv` is not tab delimited
    colnames <- read_lines(colFile)
    # Move the multiplexed sample index number to the beginning for easier
    # downstream grep matching
    colnames <- sub("^([ACGT]+)-(\\d+)$", "\\2-\\1", colnames)

    # `genes.tsv` is tab delimited
    rownames <- read_tsv(
        file = rowFile,
        col_names = c("geneID", "geneName"),
        col_types = "cc"
    )
    assert_has_rows(rownames)
    rownames <- pull(rownames, "geneID")

    counts <- readMM(matrixFile)
    assert_are_identical(length(colnames), ncol(counts))
    assert_are_identical(length(rownames), nrow(counts))
    colnames(counts) <- colnames
    rownames(counts) <- rownames
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
    list <- mcmapply(
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
