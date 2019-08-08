## FIXME Consider adding BiocParallel support here.



.import.bcbio <-  # nolint
    function(sampleDirs) {
        message("Importing counts...")
        assert(
            all(isDirectory(sampleDirs)),
            hasNames(sampleDirs)
        )

        list <- mapply(
            sampleID = names(sampleDirs),
            dir = sampleDirs,
            FUN = function(sampleID, dir) {
                counts <- .import.bcbio.mtx(dir)
                ## Prefix cell barcodes with sample identifier when we're
                ## loading counts from multiple samples.
                if (length(sampleDirs) > 1L) {
                    colnames(counts) <-
                        paste(sampleID, colnames(counts), sep = "_")
                }
                ## Ensure names are valid.
                counts <- makeDimnames(counts)
                counts
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )

        ## Remove any empty items in list, which can result from low quality
        ## samples with empty matrices in bcbio pipeline.
        list <- Filter(Negate(is.null), list)
        if (!hasLength(list)) {
            ## nocov start
            stop(paste(
                "bcbio didn't return any cells.",
                "Check your `minimum_barcode_depth` setting."
            ))
            ## nocov end
        }

        ## Bind the matrices.
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
            sampleDirs,
            paste(basename(sampleDirs), "barcodes.tsv", sep = "-")
        )
        files <- realpath(files)
        names(files) <- names(sampleDirs)
        list <- lapply(files, function(file) {
            data <- read_tsv(
                file = file,
                col_names = c("barcode", "n"),
                col_types = "ci"
            )
            x <- data[["n"]]
            names(x) <- makeNames(data[["barcode"]])
            assert(is.integer(x))
            x
        })
        names(list) <- names(sampleDirs)
        list
    }



## Import bcbio Sparse Counts Data
## Always in Matrix Market Exchange (MEX/MTX) format.
## This may be advantagenous to loading the giant combined matrix because we
## can parallelize with BiocParallel.
.import.bcbio.mtx <-  # nolint
    function(dir) {
        assert(isADirectory(dir))

        ## Require that all of the files exist, even if they are empty.
        file <- file.path(dir, paste0(basename(dir), ".mtx"))
        rownamesFile <- paste0(file, ".rownames")
        colnamesFile <- paste0(file, ".colnames")
        assert(all(isFile(c(file, rownamesFile, colnamesFile))))

        ## Attempt to load the column and rowname files first. If they're empty,
        ## skip loading the MatrixMarket file, which will error otherwise. The
        ## bcbio pipeline will output empty files for very low quality samples
        ## with no cells that pass filtering.

        ## Import Genes/transcripts (features).
        rownames <- read_lines(rownamesFile)

        ## Import cellular barcodes.
        colnames <- read_lines(colnamesFile)

        if (!length(rownames) > 0L || !length(colnames) > 0L) {
            ## nocov start
            message(paste0("Skipped", basename(dir), "."))
            return(NULL)
            ## nocov end
        }

        ## Import counts.
        counts <- readMM(file)

        assert(
            identical(length(rownames), nrow(counts)),
            identical(length(colnames), ncol(counts))
        )

        rownames(counts) <- rownames
        colnames(counts) <- colnames

        message(paste0("Imported ", basename(dir), "."))

        counts
    }
