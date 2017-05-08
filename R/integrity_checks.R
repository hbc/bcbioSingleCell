#' Perform object integrity checks
#'
#' @rdname integrity_checks
#' @author Michael Steinbaugh



#' @rdname integrity_checks
#' @description Check \code{bcbio-nextgen} run
#' @param run \code{bcbio-nextgen} run
#' @export
check_run <- function(run) {
    if (!is.list(run)) {
        stop("Run object is not a list")
    }
    if (!file.exists(run$upload_dir)) {
        stop("Could not access upload dir")
    }
    if (!length(dir(run$project_dir))) {
        stop("Could not access project_dir")
    }
    if (!length(run$sample_dirs)) {
        stop("No sample directories in run")
    }
}



#' @rdname integrity_checks
#' @description Check sparse counts matrix
#' @param sparsecounts Sparse counts matrix
#' @export
check_sparse <- function(sparsecounts) {
    if (class(sparsecounts)[1] != "dgCMatrix") {
        stop("sparse counts must be dgCMatrix")
    }
}
