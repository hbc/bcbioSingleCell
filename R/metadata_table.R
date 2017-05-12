#' Metadata table for knit report.
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#'
#' @return \code{\link[knitr]{kable}}.
#' @export
metadata_table <- function(run) {
    run$metadata %>%
        remove_rownames %>%
        mutate(file_name = NULL) %>%
        set_names(str_replace_all(colnames(.), "_", " ")) %>%
        kable(caption = "Sample metadata")
}
