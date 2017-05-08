#' Load bcbio-nextgen run
#'
#' We recommend loading the \code{bcbio-nextgen} run as a remote connection over
#' \code{sshfs}.
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running \code{bcbio_nextgen -w template}.
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes.
#' @param metadata Optional custom metadata file to import
#' @param read_counts Automatically read in the count data using
#'   \code{read_bcbio_sparsecounts()}
#'
#' @return \code{bcbio-nextgen} run object
#' @export
load_run <- function(
    upload_dir = "final",
    intgroup = "sample_name",
    organism,
    metadata,
    read_counts = TRUE) {
    if (!length(dir(upload_dir))) {
        stop("final directory failed to load")
    }
    upload_dir <- normalizePath(upload_dir)

    # project_dir
    pattern <- ".*/(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    match <- dir(upload_dir, full.names = TRUE) %>%
        str_subset(pattern) %>%
        str_match(pattern)
    project_dir <- match[1]
    run_date <- match[2]
    template <- match[3]

    message(paste(dir(upload_dir), collapse = "\n"))

    # Sample directories
    sample_dirs <- dir(upload_dir, full.names = TRUE)
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
        upload_dir = upload_dir,
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

    # Save Ensembl annotations for all genes
    run$ensembl <- ensembl_annotations(run)
    run$ensembl_version <- listMarts() %>% .[1, 2]

    # Read counts into a sparse matrix
    if (isTRUE(read_counts)) {
        run$sparse <- read_bcbio_sparsecounts(run)
    }

    check_run(run)
    return(run)
}
