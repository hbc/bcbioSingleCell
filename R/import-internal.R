#' Import bcbio counts from sample directories
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @note Updated 2021-09-03.
#' @noRd
#'
#' @inheritParams BiocParallel::bpmapply
#' @param sampleDirs `character`.
#' Sample directory paths.
.importCounts <- # nolint
    function(sampleDirs,
             BPPARAM # nolint
    ) {
        assert(
            allAreDirectories(sampleDirs),
            hasNames(sampleDirs),
            identical(attr(class(BPPARAM), "package"), "BiocParallel")
        )
        alert("Importing counts.")
        list <- bpmapply(
            sampleId = names(sampleDirs),
            dir = sampleDirs,
            FUN = function(sampleId, dir) {
                counts <- .importCountsPerSample(dir)
                ## Prefix cell barcodes with sample identifier when we're
                ## loading counts from multiple samples.
                if (length(sampleDirs) > 1L) {
                    colnames(counts) <-
                        paste(sampleId, colnames(counts), sep = "_")
                }
                ## Ensure names are valid.
                counts <- makeDimnames(counts)
                counts
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE,
            BPPARAM = BPPARAM
        )
        ## Remove any empty items in list, which can result from low quality
        ## samples with empty matrices in bcbio pipeline.
        list <- Filter(f = Negate(is.null), x = list)
        assert(
            hasLength(list),
            msg = sprintf(
                fmt = paste0(
                    "bcbio didn't return any cells.\n",
                    "Check your '%s' setting."
                ),
                "minimum_barcode_depth"
            )
        )
        ## Bind the matrices.
        do.call(cbind, list)
    }



#' Import counts per sample from sparse matrix
#'
#' Always in Matrix Market Exchange (MEX/MTX) format.
#'
#' This may be advantagenous to loading the giant combined matrix because we
#' can parallelize with BiocParallel.
#'
#' Attempt to load the column and rowname files first. If they're empty, skip
#' loading the MatrixMarket file, which will error otherwise. The bcbio pipeline
#' will output empty files for very low quality samples with no cells that pass
#' filtering.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @note Updated 2020-01-20.
#' @noRd
#'
#' @param dir `character(1)`.
#' Sample directory path.
#'
#' @return `sparseMatrix`.
.importCountsPerSample <- # nolint
    function(dir) {
        assert(isADirectory(dir))
        ## Require that all of the files exist, even if they are empty.
        file <- file.path(dir, paste0(basename(dir), ".mtx"))
        rownamesFile <- paste0(file, ".rownames")
        colnamesFile <- paste0(file, ".colnames")
        assert(allAreFiles(c(file, rownamesFile, colnamesFile)))
        ## Import Genes/transcripts (features).
        rownames <- import(rownamesFile, format = "lines")
        ## Import cellular barcodes.
        colnames <- import(colnamesFile, format = "lines")
        if (!length(rownames) > 0L || !length(colnames) > 0L) {
            ## nocov start
            alertWarning(sprintf("Skipped {.path %s}.", basename(dir)))
            return(NULL)
            ## nocov end
        }
        ## Import counts.
        counts <- import(file)
        assert(
            identical(length(rownames), nrow(counts)),
            identical(length(colnames), ncol(counts))
        )
        rownames(counts) <- rownames
        colnames(counts) <- colnames
        alert(sprintf("Imported {.path %s}.", basename(dir)))
        counts
    }



#' Import raw cellular barcode read list
#'
#' Get the number of pre-UMI disambiguated reads per cellular barcode.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @note Updated 2021-02-22.
#' @noRd
#'
#' @inheritParams BiocParallel::bplapply
#' @param sampleDirs `character`.
#' Sample directories.
#'
#' @return `list`.
#' List of integer vectors per sample containing the pre-filtered cellular
#' barcode counts (`nCount`).
.importReads <- # nolint
    function(sampleDirs,
             BPPARAM # nolint
    ) {
        assert(
            allAreDirectories(sampleDirs),
            hasNames(sampleDirs),
            identical(attr(class(BPPARAM), "package"), "BiocParallel")
        )
        alert("Importing unfiltered cellular barcode distributions.")
        files <- file.path(
            sampleDirs,
            paste(basename(sampleDirs), "barcodes.tsv", sep = "-")
        )
        files <- realpath(files)
        names(files) <- names(sampleDirs)
        list <- bplapply(
            X = files,
            FUN = function(file) {
                data <- import(
                    con = file,
                    format = "tsv",
                    colnames = c("barcode", "n")
                )
                x <- as.integer(data[["n"]])
                names(x) <- makeNames(data[["barcode"]])
                x
            },
            BPPARAM = BPPARAM
        )
        names(list) <- names(sampleDirs)
        list
    }
