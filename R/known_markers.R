#' Known markers data frame
#'
#' @param markers_file Gene markers file (CSV or Excel).
#' @param genome_build Genome build.
#'
#' @export
known_markers <- function(markers_file, genome_build) {
    annotable <- annotable(genome_build)
    read_file_by_extension(markers_file) %>%
        snake %>%
        left_join(annotable, by = "symbol") %>%
        filter(!is.na(.data[["ensgene"]]))
}
