#' Create a skeleton bcbio project
#'
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param base_dir Base directory in which to create skeleton bcbio structure
#' @param config_dir Configuration directory
#' @param upload_dir Upload directory
#' @param project_dir Project directory (YYYY-MM-DD_template), nested inside
#'   upload directory
#' @param sample_dirs Sample directories, nested inside upload directory
#' @param organism Organism to use
#'
#' @return bcbio run skeleton
#' @export
create_bcbio_skeleton <- function(
    base_dir = getwd(),
    run_dir = "skeleton",
    config_dir = "config",
    upload_dir = "final",
    project_dir = NULL,
    sample_dirs = NULL,
    organism,
    intgroup = "sample_name") {
    # run_dir
    if (!is.null(run_dir)) {
        run_dir <- file.path(base_dir, run_dir) %>% normalizePath
        dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    } else {
        stop("Run directory required")
    }

    # config_dir
    if (!is.null(config_dir)) {
        config_dir <- file.path(run_dir, config_dir)
        dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
    } else {
        stop("Configuration directory required")
    }

    # upload_dir
    if (!is.null(upload_dir)) {
        upload_dir <- file.path(run_dir, upload_dir)
        dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)
    } else {
        stop("Final upload directory required")
    }

    # project_dir
    if (!is.null(project_dir)) {
        project_dir <- file.path(upload_dir, project_dir)
        dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # sample_dirs
    if (!is.null(sample_dirs)) {
        sample_dirs <- file.path(upload_dir, sample_dirs)
        dir.create(sample_dirs, recursive = TRUE, showWarnings = FALSE)
    }

    run <- list(run_dir = run_dir,
                config_dir = config_dir,
                upload_dir = upload_dir,
                project_dir = project_dir,
                sample_dirs = sample_dirs,
                organism = organism,
                intgroup = intgroup)
    return(run)
}
