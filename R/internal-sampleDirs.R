#' Detect Sample Directories
#'
#' @rdname sampleDirs
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param uploadDir Upload directory.
#' @param pipeline Pipeline used to generate the samples.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [stop] if no complete sample directories match.
.sampleDirs <- function(uploadDir, pipeline = "bcbio") {
    if (pipeline == "bcbio") {
        sampleDirs <- list.dirs(
            uploadDir, full.names = TRUE, recursive = FALSE)
        # Remove the nested `projectDir`
        if (any(str_detect(basename(sampleDirs), projectDirPattern))) {
            sampleDirs <- sampleDirs %>%
                .[!str_detect(basename(.), projectDirPattern)]
        }
        names(sampleDirs) <- basename(sampleDirs)
    } else if (pipeline == "cellranger") {
        # Use the raw, not filtered matrices
        sampleDirs <- list.files(
            uploadDir, pattern = "^matrix\\.mtx$",
            full.names = TRUE, recursive = TRUE) %>%
            # Look for the filtered counts
            str_subset(file.path("outs", "filtered_gene_bc_matrices")) %>%
            dirname
        # Recurse through file path to get names
        # cellranger/SAMPLE/outs/filtered_gene_bc_matrices/GENOME"
        names(sampleDirs) <- sampleDirs %>%
            dirname %>%
            dirname %>%
            dirname %>%
            basename
    } else {
        stop("Unsupported pipeline")
    }

    # Return
    if (length(sampleDirs) == 0L) {
        stop("No sample directories detected")
    } else {
        message(paste(length(sampleDirs), "samples detected"))
        sampleDirs <- camel(sampleDirs)
    }
    sampleDirs
}
