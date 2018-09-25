# bcbio ========================================================================
.import.bcbio <-  # nolint
    function(sampleDirs) {
        assert_all_are_dirs(sampleDirs)
        assert_has_names(sampleDirs)
        message("Importing counts...")
        list <- mapply(
            sampleID = names(sampleDirs),
            dir = sampleDirs,
            FUN = function(sampleID, dir) {
                counts <- .import.bcbio.mtx(sampleID, dir)
            }
        )
        # Remove any empty items in list, which can result from low quality
        # samples with empty matrices in bcbio pipeline.
        list <- Filter(Negate(is.null), list)
        # Bind the matrices.
        do.call(cbind, list)
    }



#' Import Pre-Filtered Cellular Barcodes List
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param sampleDirs Sample directories.
#'
#' @return `list`.
.import.bcbio.barcodes <-  # nolint
    function(sampleDirs) {
        files <- file.path(
            normalizePath(sampleDirs, winslash = "/", mustWork = TRUE),
            paste(basename(sampleDirs), "barcodes.tsv", sep = "-")
        )
        assert_all_are_existing_files(files)
        names(files) <- names(sampleDirs)
        message("Importing unfiltered cellular barcode distributions...")
        list <- lapply(files, function(file) {
            read_tsv(
                file = file,
                col_names = c("barcode", "n"),
                col_types = "ci"
            ) %>%
                as_tibble() %>%
                mutate(
                    !!sym("barcode") :=
                        makeNames(!!sym("barcode"), unique = TRUE)
                )
        })
        names(list) <- names(sampleDirs)
        list
    }



# Import bcbio Sparse Counts Data
# Always in Matrix Market Exchange (MEX/MTX) format.
# This may be advantagenous to loading the giant combined matrix because we
# can parallelize with BiocParallel.
.import.bcbio.mtx <-  # nolint
    function(sampleID, dir) {
        assert_is_a_string(sampleID)
        assertAllAreValidNames(sampleID)
        assert_is_a_string(dir)
        assert_all_are_dirs(dir)

        # Require that all of the files exist, even if they are empty.
        matrixFile <- file.path(dir, paste0(basename(dir), ".mtx"))
        rowFile <- paste0(matrixFile, ".rownames")
        colFile <- paste0(matrixFile, ".colnames")
        assert_all_are_existing_files(c(matrixFile, rowFile, colFile))

        # Attempt to load the column and rowname files first. If they're empty,
        # skip loading the MatrixMarket file, which will error otherwise. The
        # bcbio pipeline will output empty files for very low quality samples
        # with no cells that pass filtering.

        # Import Genes/transcripts (features).
        rownames <- read_lines(rowFile)

        # Import cellular barcodes.
        colnames <- read_lines(colFile)

        if (!length(rownames) > 0L || !length(colnames) > 0L) {
            message(paste("Skipped", sampleID))
            return(NULL)
        }

        # Import counts.
        counts <- readMM(matrixFile)
        assert_is_all_of(counts, "sparseMatrix")

        assert_are_identical(length(rownames), nrow(counts))
        assert_are_identical(length(colnames), ncol(counts))

        # Ensure rows are valid.
        rownames(counts) <- rownames %>%
            makeNames(unique = TRUE)

        # Prefix cell barcodes with sample ID and make valid.
        colnames(counts) <- colnames %>%
            paste(sampleID, ., sep = "_") %>%
            makeNames(unique = TRUE)

        message(paste0("Imported ", sampleID))
        counts
    }



# Cell Ranger ==================================================================
# FIXME Need to add a simple mode, like `Seurat::Read10X()`.

.import.cellranger <-  # nolint
    function(
        sampleDirs,
        format = c("mtx", "h5"),
        filtered = TRUE
    ) {
        assert_all_are_dirs(sampleDirs)
        assert_has_names(sampleDirs)

        format <- match.arg(format)
        assert_is_a_bool(filtered)

        message(paste("Reading counts in", format, "format"))
    }



# Import Cell Ranger HDF5 Counts
# @seealso `cellrangerRkit::get_matrix_from_h5()`
.import.cellranger.h5 <-  # nolint
    function(sampleDir, filtered = TRUE) {
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

        # Want `Csparse` instead of `Tsparse` generally.
        counts <- sparseMatrix(
            i = h5[["indices"]] + 1L,
            p = h5[["indptr"]],
            x = as.numeric(h5[["data"]]),
            dims = h5[["shape"]],
            giveCsparse = TRUE
        )

        rownames <- h5[["genes"]]
        colnames <- h5[["barcodes"]]
        # Move the multiplexed sample index number to the beginning.
        colnames <- sub("^([ACGT]+)-(\\d+)$", "\\2-\\1", colnames)
        assert_are_identical(length(rownames), nrow(counts))
        assert_are_identical(length(colnames), ncol(counts))
        rownames(counts) <- rownames
        colnames(counts) <- colnames

        counts
    }



# Import Cell Ranger Sparse Counts
# Matrix Market Exchange (MEX/MTX) format.
.import.cellranger.mtx <-
    function(sampleDir, filtered = TRUE) {
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
