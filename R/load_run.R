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
    # Switch to `bcbioScDataSet` S4 class
    run <- list()

    # upload_dir
    if (!dir.exists(upload_dir)) {
        stop("Upload directory missing")
    }
    run$upload_dir <- normalizePath(upload_dir)

    # metadata_file
    if (!file.exists(metadata_file)) {
        stop("Metadata file missing")
    }
    run$metadata_file <- normalizePath(metadata_file)

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

    # sample_dirs. Subset later using metadata data frame.
    sample_dirs <- dir(run$upload_dir, full.names = TRUE) %>%
        set_names(basename(.)) %>%
        # Remove the nested `project_dir`
        .[!grepl(basename(run$project_dir), names(.))]
    message(paste(length(sample_dirs), "samples detected in run"))

    # Detect number of sample lanes
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2] %>% unique %>% length
        message(paste(lanes, "lane replicates per sample detected"))
    } else {
        lanes <- NULL
    }
    run$lanes <- lanes

    # Metadata data frame, with sample_barcode rownames
    run$metadata <- read_metadata(run)

    # Check that metadata matches the sample_dirs
    if (!all(run$metadata$sample_barcode %in% names(sample_dirs))) {
        stop("Sample name mismatch between directories and metadata")
    }

    # Subset sample dirs by matching sample barcodes in metadata
    run$sample_dirs <- sample_dirs[run$metadata$sample_barcode]
    message(paste(length(run$sample_dirs), "samples matched by metadata"))

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

    # Get run date and template from project_dir
    match <- str_match(project_dir, pattern)
    run$date <- c(bcbio = as.Date(match[2]),
                  R = Sys.Date())
    run$template <- match[3]

    # Final slots
    run$wd <- getwd()
    run$hpc <- detect_hpc()
    run$session <- sessionInfo()

    check_run(run)
    return(run)
}
