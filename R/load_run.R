#' Load bcbio-nextgen run
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param metadata_file (**Required**). Sample barcode metadata file.
#' @param organism Organism name, following Ensembl conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes.
#' @param interesting_groups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#'
#' @return [bcbioSCDataSet].
#' @export
load_run <- function(
    upload_dir = "final",
    # FIXME External metadata file is required until YAML support added
    metadata_file,
    # FIXME Can use versions file to detect genome, instead of YAML
    organism = NULL,
    interesting_groups = "sample_name") {

    # Upload directory
    if (!dir.exists(upload_dir)) {
        stop("Upload directory missing")
    }
    upload_dir <- normalizePath(upload_dir)
    pipeline <- .detect_pipeline(upload_dir)

    # Sample directories
    # Locate nested matrix.mtx files to define sample directories. A similar
    # code approach is used in [.detect_pipeline()].
    sample_dirs <- list.files(
        upload_dir, pattern = "*.mtx$",
        full.names = TRUE, recursive = TRUE) %>%
        normalizePath %>%
        dirname %>%
        set_names(basename(.))

    # bcbio-nextgen ====
    if (pipeline == "bcbio") {
        # Find most recent nested project_dir (normally only 1)
        project_dir_pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
        project_dir <- dir(upload_dir,
                           pattern = project_dir_pattern,
                           full.names = FALSE,
                           recursive = FALSE)
        if (length(project_dir) != 1) {
            stop("Uncertain about project directory location")
        }
        message(project_dir)
        match <- str_match(project_dir, project_dir_pattern)
        date <- c(bcbio = as.Date(match[[2L]]),
                  R = Sys.Date())
        template <- match[[3L]]
        project_dir <- file.path(upload_dir, project_dir)

        # Remove the nested `project_dir` from sample directories
        sample_dirs <- sample_dirs %>%
            .[!grepl(basename(project_dir), names(.))]
    }

    # 10X Chromium CellRanger ====
    if (pipeline == "cellranger") {
        if (!length(matrices)) {
            stop("No count matrices detected")
        }
        parent_dirs <- dirname(matrices)
    }

    if (!length(sample_dirs)) {
        stop("No sample directories in run")
    }
    message(paste(length(sample_dirs), "samples detected in run"))

    # Detect number of sample lanes
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2L] %>% unique %>% length
        message(paste(lanes, "lane replicates per sample detected"))
    } else {
        lanes <- 1L
    }

    # Metadata data frame, with sample_barcode rownames
    # FIXME Make sure custom metadata stops if missing
    custom_metadata <- .custom_metadata(metadata_file)

    # Check that sample directories match the metadata
    if (!all(names(sample_dirs) %in% metadata[["sample_barcode"]])) {
        stop("Sample name mismatch between directories and metadata")
    }

    # Subset sample dirs by matching sample barcodes in metadata
    sample_dirs <- sample_dirs[metadata[["sample_barcode"]]]
    message(paste(length(sample_dirs), "samples matched by metadata"))

    # Program versions
    programs <- .programs(project_dir)

    # Read counts into a sparse matrix
    message("Reading counts into sparse matrix...")
    sparse_counts <- .sparse_counts(project_dir)

    # Read cellular barcode distributions
    message("Reading cellular barcode distributions...")
    barcodes <- .barcodes(sample_dirs)

    # Generate barcode metrics
    message("Calculating barcode metrics...")
    metrics <- .metrics(sparse_counts)

    # col_data
    col_data <- custom_metadata[colnames(counts), ]
    identical(colnames(sparse_counts), rownames(colData))

    # row_data
    # FIXME use annotables...follow the bcbioRnaseq code

    counts <- NULL
    sparse <- NULL

    # Metadata ====
    metadata <- SimpleList(
        wd = getwd(),
        hpc = detect_hpc(),
        session_info = sessionInfo())

    # Package into SummarizedExperiment ====
    se <- .summarized_experiment(
        counts,
        col_data = col_data,
        row_data = row_data,
        metadata = metadata)

    bcb <- new("bcbioRnaDataSet", se)
    bcbio(bcb, "metrics") <- metrics
    bcb
}
