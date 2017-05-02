#' Metadata table for knit report
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run
#'
#' @return kable
#' @export
metadata_table <- function(run) {
    run$metadata %>%
        remove_rownames %>%
        mutate(description = NULL,
               file_name = NULL,
               reverse_complement = NULL) %>%
        set_names(str_replace_all(colnames(.), "_", " ")) %>%
        kable(caption = "Sample metadata")
}
