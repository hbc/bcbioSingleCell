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
    if (!file.exists(run$parent_dir)) {
        stop("could not access parent_dir")
    }
    if (!file.exists(run$run_dir)) {
        stop("could not access run_dir")
    }
    if (!file.exists(run$config_dir)) {
        stop("could not access config_dir")
    }
    if (!file.exists(run$final_dir)) {
        stop("could not access final_dir")
    }
    if (!length(run$sample_dirs)) {
        stop("no sample directories in run")
    }
    if (!file.exists(run$project_dir)) {
        stop("could not access project_dir")
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

##' create a skeleton bcbio project
##'
##' @param basedir base directory in which to create skeleton bcbio structure
##' @param organism organism to use
##' @return a bcbio run object
##' @author Rory Kirchner
##' @author Michael Steinbaugh
##' @export
make_bcbio_skeleton <- function(basedir=getwd(), organism="hsapiens") {
  parent_dir = file.path(basedir, "skeleton")
  run <- list(parent_dir = file.path(parent_dir),
              run_dir = file.path(parent_dir),
              config_dir = file.path(parent_dir, "config"),
              final_dir = file.path(parent_dir, "final"),
              project_dir = file.path(parent_dir, "final"),
              sample_dirs = c(NA),
              organism = organism,
              intgroup = NA)
  dir.create(run$parent_dir, showWarnings = FALSE)
  dir.create(run$config_dir, showWarnings = FALSE)
  dir.create(run$final_dir, showWarnings = FALSE)
  dir.create(run$project_dir, showWarnings = FALSE)
  return(run)
}
