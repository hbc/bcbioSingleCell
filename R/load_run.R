#' Load bcbio-nextgen run.
#'
#' We recommend loading the bcbio-nextgen run as a remote connection over
#' \code{sshfs}.
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running \code{bcbio_nextgen -w template}.
#' @param metadata_file Sample barcode metadata file.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes.
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#'
#' @return bcbio-nextgen run object.
#' @export
load_run <- function(
    upload_dir = "final",
    metadata_file = file.path("meta", "indrop_rnaseq.xlsx"),
    organism,
    intgroup = "sample_name") {
    # Check connection to upload_dir
    if (!dir.exists(upload_dir)) {
        stop("Upload directory missing")
    }
    upload_dir <- normalizePath(upload_dir)

    # Check metadata_file
    if (!file.exists(metadata_file)) {
        stop("Metadata file missing")
    }
    metadata_file <- normalizePath(metadata_file)

    # Find most recent nested project_dir (normally only 1)
    pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    project_dir <- dir(upload_dir,
                       pattern = pattern,
                       full.names = FALSE,
                       recursive = FALSE) %>%
        sort %>% rev %>% .[[1]]
    if (!length(project_dir)) {
        stop("Failed to match project directory")
    }
    message(project_dir)
    match <- str_match(project_dir, pattern)
    run_date <- match[2]
    template <- match[3]
    project_dir <- file.path(upload_dir, project_dir)
    if (!dir.exists(project_dir)) {
        stop("Project directory missing")
    }

    # Get all sample folders. We will subset this vector using the metadata
    # data frame later.
    sample_dirs <- dir(upload_dir, full.names = TRUE) %>%
        set_names(basename(.)) %>%
        # Remove the nested `project_dir`
        .[!grepl(basename(project_dir), names(.))]
    message(paste(length(sample_dirs), "samples detected in run"))

    # Detect number of sample lanes
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        message("Lane split technical replicates detected")
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2] %>% unique %>% length
    } else {
        lanes <- NULL
    }

    # Set up the run list skeleton. Switch this to an S4 class once we finish
    # updating the bcbioRnaseq package. `bcbioSCDataSet`.
    run <- list(
        upload_dir = upload_dir,
        project_dir = project_dir,
        sample_dirs = sample_dirs,
        metadata_file = metadata_file,
        organism = organism,
        intgroup = intgroup,
        lanes = lanes,
        today_date = Sys.Date(),
        run_date = as.Date(run_date),
        template = template,
        wd = getwd(),
        hpc = detect_hpc(),
        session_info = sessionInfo()
    )

    run$metadata <- read_metadata(run)
    # Subset sample dirs by matching sample barcodes in metadata
    run$sample_dirs <- run$sample_dirs[run$metadata$sample_barcode]
    if (length(run$sample_dirs)) {
        message(paste(length(sample_dirs), "samples matched by metadata"))
    }

    # Save Ensembl annotations for all genes
    run$ensembl <- ensembl_annotations(run)
    run$ensembl_version <- listMarts() %>% .[1, 2]

    # Program versions
    message("Reading program versions...")
    run$programs <- file.path(project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))

    # Read counts into a sparse matrix
    message("Reading counts into sparse matrix...")
    run$counts <- read_bcbio_counts(run)

    # Generate barcode metrics
    message("Calculating barcode metrics...")
    run$metrics <- barcode_metrics(run)

    check_run(run)
    return(run)
}
