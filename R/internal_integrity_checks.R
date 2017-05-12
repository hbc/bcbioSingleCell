#' Perform object integrity checks.
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
    if (!dir.exists(run$upload_dir)) {
        stop("Could not access upload dir")
    }
    if (!dir.exists(run$project_dir)) {
        stop("Could not access project_dir")
    }
    if (!all(dir.exists(run$sample_dirs))) {
        stop("No sample directories in run")
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
