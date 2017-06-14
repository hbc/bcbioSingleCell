#' Check object integrity
#'
#' @rdname integrity
#' @keywords internal
#' @param run [bcbioSCDataSet].
#' @export
check_run <- function(run) {
    if (!is.list(run)) {
        stop("Run object is not a list")
    }
}



#' @rdname integrity
#' @param sparse Sparse counts matrix.
#' @export
check_sparse <- function(sparse) {
    if (class(sparse)[1] != "dgCMatrix") {
        stop("Sparse counts must be dgCMatrix")
    }
}



#' Create a skeleton bcbio project
#'
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param base_dir Base directory in which to create skeleton bcbio structure.
#' @param config_dir Configuration directory.
#' @param upload_dir Upload directory.
#' @param project_dir Project directory (YYYY-MM-DD_template), nested inside
#'   upload directory.
#' @param sample_dirs Sample directories, nested inside upload directory.
#' @param organism Organism to use.
#'
#' @return bcbio run skeleton.
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

    list(run_dir = run_dir,
         config_dir = config_dir,
         upload_dir = upload_dir,
         project_dir = project_dir,
         sample_dirs = sample_dirs,
         organism = organism,
         intgroup = intgroup)
}



load_run_as_list <- function(
    upload_dir,
    metadata_file,
    organism,
    intgroup) {
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
    run$programs <- read_bcbio_programs(run)

    # Read counts into a sparse matrix
    message("Reading counts into sparse matrix...")
    run$sparse <- read_bcbio_counts(run)

    # Read cellular barcode distributions
    message("Reading cellular barcode distributions...")
    run$barcodes <- read_bcbio_barcodes(run)

    # Generate barcode metrics
    message("Calculating barcode metrics...")
    run$metrics <- metrics(run)

    # Final slots
    run$wd <- getwd()
    run$hpc <- detect_hpc()
    run$session <- sessionInfo()

    run
}
