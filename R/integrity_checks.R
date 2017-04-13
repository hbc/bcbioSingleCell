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
        stop("bcbio run object not found")
    }
    if (!length(dir(run$parent_dir))) {
        stop("could not access parent_dir")
    }
    if (!length(dir(run$run_dir))) {
        stop("could not access run_dir")
    }
    if (!length(dir(run$config_dir))) {
        stop("could not access config_dir")
    }
    if (!length(dir(run$final_dir))) {
        stop("could not access final_dir")
    }
    if (!length(run$sample_dirs)) {
        stop("no sample directories in run")
    }
    if (!length(dir(run$project_dir))) {
        stop("could not access project_dir")
    }
}
