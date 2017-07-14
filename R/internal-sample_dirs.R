#' Detect sample directories
#'
#' @rdname sample_dirs
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Upload directory.
#' @param pipeline Pipeline used to generate the samples.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [stop] if no complete sample directories match.
.sample_dirs <- function(upload_dir, pipeline = "bcbio") {
    if (pipeline == "bcbio") {
        sample_dirs <- list.dirs(
            upload_dir, full.names = TRUE, recursive = FALSE)
        # Remove the nested `project_dir`
        if (any(str_detect(basename(sample_dirs), project_dir_pattern))) {
            sample_dirs <- sample_dirs %>%
                .[!str_detect(basename(.), project_dir_pattern)]
        }
    } else if (pipeline == "cellranger") {
        # Use the raw, not filtered matrices
        sample_dirs <- list.files(
            upload_dir, pattern = "^matrix\\.mtx$",
            full.names = TRUE, recursive = TRUE) %>%
            # Look for the filtered counts
            str_subset(file.path("outs", "filtered_gene_bc_matrices")) %>%
            dirname
    } else {
        stop("Unsupported pipeline")
    }

    # Return
    if (length(sample_dirs) == 0L) {
        stop("No sample directories detected")
    } else {
        # Generate names from file paths
        names <- sample_dirs %>%
            # Strip upload_dir
            str_replace(str_c(upload_dir, "/"), "") %>%
            # Strip subfolders
            #str_replace("/.*", "") %>%
            # Make valid names
            make.names
        sample_dirs <- normalizePath(sample_dirs)
        names(sample_dirs) <- names
        message(paste(length(sample_dirs), "samples detected"))
    }
    sample_dirs
}
