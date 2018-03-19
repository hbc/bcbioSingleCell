#' Subset Per Sample
#'
#' @name subsetPerSample
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param minCells Minimum number of cells required per sample.
#' @param envir Environment where to assign the subsets.
#' @param dir Output directory.
#'
#' @return Named vector of saved `bcbioSingleCell` subset file paths.
#'
#' @examples
#' # bcbioRNASeq ====
#' # This will save the subsets per sample to disk
#' subsetPerSample(bcb_small, dir = "subsetPerSample")
#' list.files("subsetPerSample")
#'
#' # Clean up
#' unlink("subsetPerSample", recursive = TRUE)
NULL



# Methods ======================================================================
#' @rdname subsetPerSample
#' @importFrom basejump assignAndSaveData initializeDirectory
#' @importFrom dplyr pull
#' @importFrom pbapply pblapply
#' @export
setMethod(
    "subsetPerSample",
    signature("bcbioSingleCell"),
    function(
        object,
        minCells = 200L,
        envir = parent.frame(),
        dir = "."
    ) {
        assertIsAnImplicitInteger(minCells)
        assert_all_are_positive(minCells)
        assert_is_environment(envir)
        dir <- initializeDirectory(dir)
        sampleIDs <- sampleData(object) %>%
            pull("sampleID") %>%
            as.character()
        files <- mapply(
            FUN = function(object, sampleID, minCells, envir, dir) {
                subset <- selectSamples(object, sampleID = sampleID)
                # Skip if subset doesn't have enough cells
                if (ncol(subset) < minCells) {
                    warn(paste(sampleID, "didn't pass minimum cell cutoff"))
                    return(NULL)
                }
                assignAndSaveData(
                    name = sampleID,
                    object = subset,
                    dir = dir
                )
            },
            sampleID = sampleIDs,
            MoreArgs = list(
                object = object,
                minCells = minCells,
                envir = envir,
                dir = dir
            ),
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )
        files <- Filter(Negate(is.null), files)
        names <- names(files)
        files <- unlist(files)
        files <- normalizePath(files, winslash = "/", mustWork = TRUE)
        names(files) <- names
        files
    }
)
