#' Load bcbio-nextgen run
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param sample_metadata_file (**Required**). Sample barcode metadata file.
#' @param well_metadata_file (*Optional*). Well identifier metadata file.
#' @param interesting_groups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#'
#' @return [bcbioSCDataSet].
#' @export
load_run <- function(
    upload_dir = "final",
    sample_metadata_file,
    well_metadata_file = NULL,
    interesting_groups = "sample_name") {


    # Upload directory ====
    if (!dir.exists(upload_dir)) {
        stop("Upload directory missing")
    }
    upload_dir <- normalizePath(upload_dir)


    # Project directory ====
    project_dir <- dir(upload_dir,
                       pattern = project_dir_pattern,
                       full.names = FALSE,
                       recursive = FALSE)
    if (length(project_dir) != 1L) {
        stop("Uncertain about project directory location")
    }
    message(project_dir)
    match <- str_match(project_dir, project_dir_pattern)
    run_date <- match[[2L]] %>% as.Date
    template <- match[[3L]]
    project_dir <- file.path(upload_dir, project_dir)


    # Log files ====
    message("Reading log files")
    bcbio_nextgen <- read_lines(
        file.path(project_dir, "bcbio-nextgen.log"))
    bcbio_nextgen_commands <- read_lines(
        file.path(project_dir, "bcbio-nextgen-commands.log"))


    # Data versions and programs ====
    data_versions <- .data_versions(project_dir)
    programs <- .programs(project_dir)
    genome_build <- data_versions %>%
        as.data.frame %>%
        filter(.data[["resource"]] == "transcripts") %>%
        pull("genome")


    # Molecular barcode (UMI) type ====
    message("Detecting UMI type")
    if (any(str_detect(bcbio_nextgen_commands,
                       "work/umis/harvard-indrop-v3.json"))) {
        umi_type <- "indrop_v3"
    }
    message(umi_type)


    # Sample directories ====
    sample_dirs <- .sample_dirs(upload_dir)


    # Sample metadata ====
    sample_metadata_file <- normalizePath(sample_metadata_file)
    sample_metadata <- .sample_metadata(sample_metadata_file, sample_dirs)


    # Well metadata ====
    if (!is.null(well_metadata_file)) {
        well_metadata_file <- normalizePath(well_metadata_file)
    }
    well_metadata <- .read_file(well_metadata_file)
    # FIXME Need a working example from Rory
    # Ensure `sample_id` is factor, by = "well_id"


    # Prepare samples based on metadata ====
    # Check for sample directory match based on metadata
    if (!all(sample_metadata[["sample_barcode"]] %in% names(sample_dirs))) {
        stop("Sample directory names don't match the metadata file")
    }

    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (!identical(sample_metadata[["sample_barcode"]], names(sample_dirs)) &
        length(sample_metadata[["sample_barcode"]]) < length(sample_dirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        all_samples <- FALSE
        sample_dirs <- sample_dirs %>%
            .[names(sample_dirs) %in% sample_metadata[["sample_barcode"]]]
        message(paste(length(sample_dirs), "samples matched by metadata"))
    } else {
        all_samples <- TRUE
    }

    # Finally, check that sample directories match the metadata
    if (!identical(sample_metadata[["sample_barcode"]], names(sample_dirs))) {
        stop("Sample name mismatch between directories and metadata")
    }


    # Sparse counts ====
    # The pipeline outputs transcript-level counts. Let's store these inside the
    # bcbioSCDataSet but outside of the SummarizedExperiment, where we will
    # instead slot the gene-level counts.
    tx_sparse_list <- pblapply(seq_along(sample_dirs), function(a) {
        .sparse_counts(sample_dirs[a])
    }
    ) %>% set_names(names(sample_dirs))
    message("Combining counts into a single sparse matrix")
    tx_sparse_counts <- do.call(cBind, tx_sparse_list)
    rm(tx_sparse_list)
    sparse_counts <- .sparse_counts_tx2gene(tx_sparse_counts, genome_build)


    # Cellular barcodes ====
    cellular_barcodes <- .cellular_barcodes(sample_dirs)


    # Row data ====
    annotable <- .annotable(genome_build)


    # Column data ====
    metrics <- .metrics(sparse_counts, annotable) %>%
        as.data.frame %>%
        rownames_to_column %>%
        left_join(.bind_cellular_barcodes(cellular_barcodes),
                  by = "rowname") %>%
        tidy_select("reads", everything()) %>%
        column_to_rownames %>%
        as.matrix


    # Metadata ====
    metadata <- SimpleList(
        upload_dir = upload_dir,
        project_dir = project_dir,
        sample_dirs = sample_dirs,
        interesting_groups = interesting_groups,
        sample_metadata_file = sample_metadata_file,
        sample_metadata = sample_metadata,
        well_metadata_file = well_metadata_file,
        well_metdata = well_metadata,
        template = template,
        data_versions = data_versions,
        programs = programs,
        bcbio_nextgen = bcbio_nextgen,
        bcbio_nextgen_commands = bcbio_nextgen_commands,
        run_date = run_date,
        load_date = Sys.Date(),
        umi_type = umi_type,
        all_samples = all_samples,
        wd = getwd(),
        hpc = detect_hpc(),
        session_info = sessionInfo())


    # Package into SummarizedExperiment ====
    se <- .summarized_experiment(
        sparse_counts,
        col_data = metrics,
        row_data = annotable,
        metadata = metadata)
    bcb <- new("bcbioSCDataSet", se)
    bcbio(bcb, "tx_sparse_counts") <- tx_sparse_counts
    bcbio(bcb, "cellular_barcodes") <- cellular_barcodes
    bcb
}
