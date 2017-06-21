#' Read inDrop i5 index counts
#'
#' @rdname indrop
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param log_file Log file produced by `indrop_i5_index_counts.sh` bash script.
#'
#' @return i5 index (barcode) counts tibble.
.indrop_i5_index_counts <- function(
    log_file = file.path("logs", "indrop_i5_index_counts.log")) {
    if (file.exists(log_file)) {
        read_table(log_file, col_names = c("counts", "reverse_complement"))
    }
}
