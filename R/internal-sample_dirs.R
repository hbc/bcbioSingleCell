# Locate nested matrix.mtx files to define sample directories. A similar
# code approach is used in [.detect_pipeline()].
.sample_dirs <- function(parent_dir) {
    sample_dirs <- list.files(
        parent_dir, pattern = "*.mtx$",
        full.names = TRUE, recursive = TRUE) %>%
        normalizePath %>%
        dirname %>%
        set_names(basename(.))
    if (!length(sample_dirs)) {
        stop("No sample directories found")
    }
    sample_dirs
}
