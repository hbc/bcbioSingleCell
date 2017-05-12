#' Read inDrop i5 index counts.
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param log_file Log file produced by \href{https://github.com/steinbaugh/seqcloud/blob/dev/bash/indrop_i5_index_counts.sh}{indrop_i5_index_counts.sh} bash script.
#'
#' @return i5 index (barcode) counts tibble.
#' @export
indrop_i5_index_counts <- function(
    log_file = file.path("logs", "indrop_i5_index_counts.log")) {
    if (file.exists(log_file)) {
        i5_index <- read_table(log_file,
                               col_names = c("counts", "reverse_complement"))
        return(i5_index)
    }
}
