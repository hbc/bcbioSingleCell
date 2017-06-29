#' Detect sample directories
#'
#' @rdname sample_dirs
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Upload directory.
#' @param pattern Match sample directories by the presence of a nested file
#'   pattern.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [stop] if no sample directories are found.
.sample_dirs <- function(upload_dir, pattern = "*\\.mtx$") {
    if (!is.null(pattern)) {
        sample_dirs <- list.files(
            upload_dir, pattern = pattern,
            full.names = TRUE, recursive = TRUE) %>%
            dirname
    } else {
        sample_dirs <- list.dirs(
            upload_dir, full.names = TRUE, recursive = FALSE)
    }
    sample_dirs <- sample_dirs %>%
        normalizePath %>%
        set_names(basename(.))

    # Remove the nested `project_dir` from a bcbio run
    if (any(str_detect(basename(sample_dirs), project_dir_pattern))) {
        sample_dirs <- sample_dirs %>%
            .[!str_detect(names(.), project_dir_pattern)]
    }
    message(paste(length(sample_dirs), "samples detected"))

    if (!length(sample_dirs)) {
        stop("No sample directories found")
    }
    sample_dirs
}
