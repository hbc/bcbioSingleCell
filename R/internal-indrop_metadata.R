# FIXME Rory applied some updates to the MASTER branch...merge here
#' Prepare inDrop metadata
#'
#' Appends reverse complement sequences by matching against the `sequence`
#' column.
#'
#' @rdname indrop_metadata
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param metadata Metadata data frame.
#'
#' @return Metadata data frame.
.indrop_metadata <- function(metadata) {
    metadata %>%
        mutate(revcomp = vapply(.data[["sequence"]], revcomp, "character"),
               sample_barcode = paste(.data[["description"]],
                                      .data[["revcomp"]],
                                      sep = "-")) %>%
        as.data.frame %>%
        set_rownames(.$sample_barcode)
}
