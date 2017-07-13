#' Read known markers
#'
#' @param markers_file Gene markers file (CSV or Excel).
#' @param genome_build Genome build.
#' @param show Show [kable].
#'
#' @return [tibble].
#' @export
read_known_markers <- function(markers_file, genome_build, show = TRUE) {
    annotable <- annotable(genome_build)
    markers <- read_file_by_extension(markers_file) %>%
        snake %>%
        left_join(annotable, by = "symbol") %>%
        filter(!is.na(.data[["ensgene"]])) %>%
        tidy_select(c("cell_type", "symbol")) %>%
        distinct
    if (isTRUE(show)) {
        kable(markers, caption = "Known markers") %>% show
    }
    markers
}
