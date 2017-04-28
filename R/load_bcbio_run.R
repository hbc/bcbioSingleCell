#' Load bcbio run
#'
#' We recommend loading the \code{bcbio-nextgen} run as a remote connection over
#' \code{sshfs}. This requires setting \code{parent_dir} in the function.
#'
#' @author Michael Steinbaugh
#'
#' @param final_dir Path to final output directory. This path is set when
#'   running \code{bcbio_nextgen -w template}.
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes.
#' @param metadata Optional custom metadata file to import
#'
#' @return \code{bcbio-nextgen} run object
#' @export
load_bcbio_run <- function(
    final_dir = "final",
    intgroup = "sample_name",
    organism,
    metadata) {
    if (!length(dir(final_dir))) {
        stop("final directory failed to load")
    }
    final_dir <- normalizePath(final_dir)

    # project_dir
    pattern <- ".*/(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    match <- dir(final_dir, full.names = TRUE) %>%
        str_subset(pattern) %>%
        str_match(pattern)
    project_dir <- match[1]
    run_date <- match[2]
    template <- match[3]

    message(paste(dir(final_dir), collapse = "\n"))

    # Sample directories
    sample_dirs <- dir(final_dir, full.names = TRUE)
    names(sample_dirs) <- basename(sample_dirs)
    # Remove the nested `project_dir`
    sample_dirs <- sample_dirs[!grepl(basename(project_dir),
                                      names(sample_dirs))]

    # Detect number of sample lanes
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2] %>% unique %>% length
    } else {
        lanes <- NULL
    }

    # Custom metadata
    if (!is.null(metadata)) {
        metadata <- read_metadata(metadata, lanes = lanes)
    }

    # Program versions
    programs <- file.path(project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))

    run <- list(
        final_dir = final_dir,
        project_dir = project_dir,
        run_date = as.Date(run_date),
        today_date = Sys.Date(),
        wd = getwd(),
        hpc = detect_hpc(),
        template = template,
        sample_dirs = sample_dirs,
        organism = organism,
        intgroup = intgroup,
        metadata = metadata,
        lanes = lanes,
        programs = programs,
        session_info = sessionInfo())

    check_run(run)
    return(run)
}
