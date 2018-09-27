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
                counts <- .import.bcbio.mtx(dir)
                # Prefix cell barcodes with sample identifier when we're loading
                # counts from multiple samples.
                if (length(sampleDirs) > 1L) {
                    colnames(counts) <-
                        paste(sampleID, colnames(counts), sep = "_")
                }
                # Ensure names are valid.
                counts <- makeDimnames(counts)
                counts
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )

        # Remove any empty items in list, which can result from low quality
        # samples with empty matrices in bcbio pipeline.
        list <- Filter(Negate(is.null), list)
        if (!has_length(list)) {
            stop(paste(
                "bcbio didn't return any cells.",
                "Check your `minimum_barcode_depth` setting."
            ))
        }

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
#' @return `list`. List of integer vectors per sample containing the
#'   pre-filtered cellular barcode counts (`nCount`).
.import.bcbio.barcodes <-  # nolint
    function(sampleDirs) {
        message("Importing unfiltered cellular barcode distributions...")
        files <- file.path(
            normalizePath(sampleDirs, winslash = "/", mustWork = TRUE),
            paste(basename(sampleDirs), "barcodes.tsv", sep = "-")
        )
        assert_all_are_existing_files(files)
        names(files) <- names(sampleDirs)
        list <- lapply(files, function(file) {
            data <- read_tsv(
                file = file,
                col_names = c("barcode", "n"),
                col_types = "ci"
            )
            x <- data[["n"]]
            names(x) <- makeNames(data[["barcode"]])
            assert_is_integer(x)
            x
        })
        names(list) <- names(sampleDirs)
        list
    }



# Import bcbio Sparse Counts Data
# Always in Matrix Market Exchange (MEX/MTX) format.
# This may be advantagenous to loading the giant combined matrix because we
# can parallelize with BiocParallel.
.import.bcbio.mtx <-  # nolint
    function(dir) {
        assert_is_a_string(dir)
        assert_all_are_dirs(dir)

        # Require that all of the files exist, even if they are empty.
        file <- file.path(dir, paste0(basename(dir), ".mtx"))
        rownamesFile <- paste0(file, ".rownames")
        colnamesFile <- paste0(file, ".colnames")
        assert_all_are_existing_files(c(file, rownamesFile, colnamesFile))

        # Attempt to load the column and rowname files first. If they're empty,
        # skip loading the MatrixMarket file, which will error otherwise. The
        # bcbio pipeline will output empty files for very low quality samples
        # with no cells that pass filtering.

        # Import Genes/transcripts (features).
        rownames <- read_lines(rownamesFile)

        # Import cellular barcodes.
        colnames <- read_lines(colnamesFile)

        if (!length(rownames) > 0L || !length(colnames) > 0L) {
            message(paste("Skipped", basename(dir)))
            return(NULL)
        }

        # Import counts.
        counts <- readMM(file)

        assert_are_identical(length(rownames), nrow(counts))
        assert_are_identical(length(colnames), ncol(counts))

        rownames(counts) <- rownames
        colnames(counts) <- colnames

        message(paste0("Imported ", basename(dir), "."))

        counts
    }



# Cell Ranger ==================================================================
.import.cellranger <-  # nolint
    function(sampleFiles) {
        assert_all_are_existing_files(sampleFiles)
        assert_has_names(sampleFiles)

        message("Importing counts...")

        if (all(grepl("\\.mtx$", sampleFiles))) {
            fun <- .import.cellranger.mtx
        } else if (all(grepl("\\.h5$", sampleFiles))) {
            fun <- .import.cellranger.h5
        }

        list <- mapply(
            sampleID = names(sampleFiles),
            file = sampleFiles,
            FUN = function(sampleID, file) {
                counts <- fun(file)
                # Prefix cell barcodes with sample identifier when we're loading
                # counts from multiple samples.
                if (length(sampleFiles) > 1L) {
                    colnames(counts) <-
                        paste(sampleID, colnames(counts), sep = "_")
                }
                # Ensure names are valid.
                counts <- makeDimnames(counts)
                counts
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )

        # Bind the matrices.
        do.call(cbind, list)
    }



# Import Cell Ranger HDF5 Counts
# @seealso `cellrangerRkit::get_matrix_from_h5()`
.import.cellranger.h5 <-  # nolint
    function(file) {
        assert_is_a_string(file)
        assert_all_are_existing_files(file)

        # Get the genome build, which we need to pass as "name".
        genomeBuild <- names(h5dump(file, load = FALSE))
        assert_is_a_string(genomeBuild)

        # Use genome build name (e.g. "/GRCh38/data").
        h5 <- h5read(file = file, name = genomeBuild)

        # Want `Csparse` not `Tsparse` matrix.
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
    function(file) {
        assert_is_a_string(file)
        assert_all_are_existing_files(file)

        # Locate required sidecar files.
        barcodesFile <- file.path(dirname(file), "barcodes.tsv")
        genesFile <- file.path(dirname(file), "genes.tsv")
        assert_all_are_existing_files(c(barcodesFile, genesFile))

        # `genes.tsv` is tab delimited.
        rownames <- read_tsv(
            file = genesFile,
            col_names = c("geneID", "geneName"),
            col_types = "cc"
        ) %>%
            pull("geneID")

        # `barcodes.tsv` is not tab delimited.
        # Move the multiplexed sample index number to the beginning for logical
        # sorting and consistency with bcbio approach.
        colnames <- read_lines(barcodesFile) %>%
            sub("^([ACGT]+)-(.+)$", "\\2-\\1", .)

        counts <- readMM(file)

        assert_are_identical(length(rownames), nrow(counts))
        assert_are_identical(length(colnames), ncol(counts))

        rownames(counts) <- rownames
        colnames(counts) <- colnames

        counts
    }
