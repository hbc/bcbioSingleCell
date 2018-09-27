# Get the MEX sparse matrix or HDF5 files inside the upload directory.
.sampleFiles.cellranger <-  # nolint
    function(
        uploadDir,
        format = c("mtx", "hdf5"),
        filtered = TRUE
    ) {
        assert_is_a_string(uploadDir)
        assert_all_are_dirs(uploadDir)
        format <- match.arg(format)
        assert_is_a_bool(filtered)

        # Look for simple upload structure.
        if ("matrix.mtx" %in% list.files(uploadDir)) {
            message("Simple mode enabled.")
            file <- file.path(uploadDir, "matrix.mtx")
            names(file) <- basename(uploadDir)
            return(file)
        }

        dirs <- list.dirs(uploadDir, recursive = FALSE)
        assert_is_non_empty(dirs)

        # Sample directories must contain `outs/` subdirectory.
        hasOuts <- vapply(
            X = dirs,
            FUN = function(dir) {
                dir.exists(file.path(dir, "outs"))
            },
            FUN.VALUE = logical(1L)
        )
        dirs <- dirs[hasOuts]
        assert_is_non_empty(dirs)

        if (isTRUE(filtered)) {
            prefix <- "filtered"
        } else {
            prefix <- "raw"
        }

        if (format == "mtx") {
            # We need to parse the subdirs, which contain the sample matrix
            # nested in a genome build directory (e.g. GRCh38, hg19).
            subdirs <- file.path(
                dirs,
                "outs",
                paste0(prefix, "_gene_bc_matrices")
            )
            assert_all_are_dirs(subdirs)
            # Get the genome build from the first sample directory.
            genomeBuild <- list.dirs(
                path = subdirs[[1L]],
                full.names = FALSE,
                recursive = FALSE
            )
            assert_is_a_string(genomeBuild)
            files <- file.path(subdirs, genomeBuild, "matrix.mtx")
        } else if (format == "hdf5") {
            files <- file.path(
                dirs,
                "outs",
                paste0(prefix, "_gene_bc_matrices_h5.h5")
            )
        }

        assert_all_are_existing_files(files)
        names(files) <- makeNames(basename(dirs))

        files
    }
