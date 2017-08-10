#' Detect Sample Directories
#'
#' @rdname internal-sampleDirs
#' @keywords internal
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
        split <- basename(sampleDirs) %>%
            str_split("-")
        sampleNames <-
            sapply(seq_along(split), function(a) {
            paste(camel(split[[a]][[1L]]),
                  split[[a]][[2L]],
                  sep = "_")
        })
        names(sampleDirs) <- sampleNames
    } else if (pipeline == "cellranger") {
        # Use the raw, not filtered matrices
        sampleDirs <- list.files(
            uploadDir, pattern = "^matrix\\.mtx$",
            full.names = TRUE, recursive = TRUE) %>%
            # Look for the filtered counts
            str_subset(file.path("filtered_gene_bc_matrices")) %>%
            dirname
        # Recurse through file path to get names
        # `cellranger/SAMPLE_ID/outs/filtered_gene_bc_matrices/GENOME`
        names(sampleDirs) <- sampleDirs %>%
            dirname %>%
            dirname %>%
            dirname %>%
            basename %>%
            camel
    } else {
        stop("Unsupported pipeline")
    }

    # Return
    if (length(sampleDirs) == 0L) {
        stop("No sample directories detected")
    }
    message(paste(length(sampleDirs), "samples detected"))
    sampleDirs
}
