#' Create a skeleton bcbio project
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param base_dir Base directory in which to create skeleton bcbio structure
#' @param organism Organism to use
#'
#' @return bcbio run skeleton
#' @export
create_bcbio_skeleton <- function(base_dir = getwd(), organism) {
    parent_dir = file.path(base_dir, "skeleton")
    run <- list(parent_dir = file.path(parent_dir),
                run_dir = file.path(parent_dir),
                config_dir = file.path(parent_dir, "config"),
                final_dir = file.path(parent_dir, "final"),
                project_dir = file.path(parent_dir, "final"),
                sample_dirs = c(NA),
                organism = organism,
                intgroup = NA)
    dir.create(run$parent_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(run$config_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(run$final_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(run$project_dir, recursive = TRUE, showWarnings = FALSE)
    return(run)
}
