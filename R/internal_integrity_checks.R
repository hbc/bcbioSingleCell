#' Perform object integrity checks
#'
#' @rdname integrity_checks
#' @author Michael Steinbaugh



#' @rdname integrity_checks
#' @keywords internal
#' @description Check bcbio-nextgen run.
#' @param run bcbio-nextgen run.
#' @export
check_run <- function(run) {
    if (!is.list(run)) {
        stop("Run object is not a list")
    }
}



#' @rdname integrity_checks
#' @keywords internal
#' @description Check sparse counts matrix.
#' @param sparse Sparse counts matrix.
#' @export
check_sparse <- function(sparse) {
    if (class(sparse)[1] != "dgCMatrix") {
        stop("sparse counts must be dgCMatrix")
    }
}
