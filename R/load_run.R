#' Load bcbio-nextgen run
#'
#' We recommend loading the bcbio-nextgen run as a remote connection over
#' \code{sshfs}.
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running \code{bcbio_nextgen -w template}.
#' @param metadata_file Sample barcode metadata file
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes.
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#' @param read_counts Automatically read in the count data using
#'   \code{read_bcbio_counts()}
#'
#' @return bcbio-nextgen run object
#' @export
load_run <- function(
    upload_dir = "final",
    metadata_file = file.path("meta", "indrop_rnaseq.xlsx"),
    organism,
    intgroup = "sample_name",
    read_counts = TRUE) {
    # Check connection to upload_dir
    if (!length(dir(upload_dir))) {
        stop("Final upload directory failed to load")
    }
    upload_dir <- normalizePath(upload_dir)

    # Find most recent nested project_dir (normally only 1)
    pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    project_dir <- dir(upload_dir,
                       pattern = pattern,
                       full.names = FALSE,
                       recursive = FALSE) %>%
        sort %>% rev %>% .[[1]]
    if (!length(project_dir)) {
        stop("Project directory not found")
    }
    message(project_dir)
    match <- str_match(project_dir, pattern)
    run_date <- match[2]
    template <- match[3]
    project_dir <- file.path(upload_dir, project_dir)

    # Program versions
    message("Reading program versions...")
    programs <- file.path(project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))

    # Obtain list of all sample folders
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

    metadata <- read_metadata(metadata_file, lanes = lanes)
    sample_dirs <- sample_dirs[metadata$sample_barcode]
    names(sample_dirs) %>% toString %>% message

    run <- list(
        upload_dir = upload_dir,
        project_dir = project_dir,
        run_date = as.Date(run_date),
        today_date = Sys.Date(),
        template = template,
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
        message("Reading counts into sparse matrix...")
        run$counts <- read_bcbio_counts(run)
        message("Reading barcode metrics...")
        run$metrics <- barcode_metrics(run, run$counts)
    }

    check_run(run)
    return(run)
}
