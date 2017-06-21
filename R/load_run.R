#' Load bcbio-nextgen run
#'
#' We recommend loading the bcbio-nextgen run as a remote connection over
#' \code{sshfs}.
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running \code{bcbio_nextgen -w template}.
#' @param sample_metadata_file Sample metadata file.
#' @param well_metadata_file Well metadata file.
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
    sample_metadata_file = NULL,
    well_metadata_file = NULL,
    organism,
    intgroup = "sample_name") {
    # Switch to `bcbioScDataSet` S4 class
    run <- list()

    # upload_dir
    if (!dir.exists(upload_dir)) {
        stop("Upload directory missing")
    }
    run$upload_dir <- normalizePath(upload_dir)

    # set sample data file if available
    run$sample_metadata_file <- sample_metadata_file
    if(!is.null(run$sample_metadata_file)) {
      if(!file.exists(run$sample_metadata_file)) {
        stop("Sample metadata file is specified, but missing.")
      }
      else {
        run$sample_metadata_file = normalizePath(run$sample_metadata_file)
      }
    }

    # set well data file if available
    run$well_metadata_file <- well_metadata_file
    if(!is.null(run$well_metadata_file)) {
      if(!file.exists(run$well_metadata_file)) {
        stop("Well metadata file is specified, but missing.")
      }
      else {
        run$sample_metadata_file = normalizePath(run$sample_metadata_file)
      }
    }

    # organism
    run$organism <- organism

    # intgroup
    run$intgroup <- intgroup

    # Find most recent nested project_dir (normally only 1)
    project_pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    project_dir <- dir(run$upload_dir,
                       pattern = project_pattern,
                       full.names = FALSE,
                       recursive = FALSE) %>%
        sort %>% rev %>% .[[1]]
    if (!length(project_dir)) {
        stop("Failed to match project directory")
    }
    message(project_dir)
    run$project_dir <- file.path(run$upload_dir, project_dir)

    # Get run date and template from project_dir
    match <- str_match(project_dir, project_pattern)
    run$date <- c(bcbio = as.Date(match[2]),
                  R = Sys.Date())
    run$template <- match[3]

    # sample_dirs. Subset later using metadata data frame.
    sample_dirs <- dir(run$upload_dir, full.names = TRUE) %>%
        set_names(basename(.)) %>%
        # Remove the nested `project_dir`
        .[!grepl(basename(run$project_dir), names(.))]
    if (!length(sample_dirs)) {
        stop("No sample directories in run")
    }
    message(paste(length(sample_dirs), "samples detected in run"))
    run$sample_dirs = sample_dirs

    # Save Ensembl annotations for all genes
    run$ensembl <- ensembl_annotations(run)
    run$ensembl_version <- listMarts() %>% .[1, 2]

    # Program versions
    message("Reading program versions...")
    run$programs <- read_bcbio_programs(run)

    # Read counts into a sparse matrix
    message("Reading counts into sparse matrix...")
    run$counts <- read_bcbio_counts(run)

    # Reading metadata about samples and wells
    run$metadata <- read_metadata(run)

    # Read cellular barcode distributions
    message("Reading cellular barcode distributions...")
    run$barcodes <- read_bcbio_barcodes(run)

    # Generate barcode metrics
    message("Calculating barcode metrics...")
    run$metrics <- barcode_metrics(run)

    # Final slots
    run$wd <- getwd()
    run$hpc <- detect_hpc()
    run$session <- sessionInfo()

    check_run(run)
    create_project_dirs()
    run
}
