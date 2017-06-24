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
#' @param interesting_groups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#'
#' @return [bcbioSCDataSet].
#' @export
load_run <- function(
    upload_dir = "final",
    metadata_file,
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
    if (length(project_dir) != 1) {
        stop("Uncertain about project directory location")
    }
    message(project_dir)
    match <- str_match(project_dir, project_dir_pattern)
    run_date <- match[[2L]] %>% as.Date
    template <- match[[3L]]
    project_dir <- file.path(upload_dir, project_dir)


    # Sample directories ====
    sample_dirs <- .sample_dirs(upload_dir)


    # Data versions and programs ====
    data_versions <- .data_versions(project_dir)
    programs <- .programs(project_dir)
    genome_build <- data_versions %>%
        as.data.frame %>%
        filter(.data[["resource"]] == "transcripts") %>%
        pull("genome")


    # Log files ====
    message("Reading log files")
    bcbio_nextgen <- read_lines(
        file.path(project_dir, "bcbio-nextgen.log"))
    bcbio_nextgen_commands <- read_lines(
        file.path(project_dir, "bcbio-nextgen-commands.log"))


    # Molecular barcode (UMI) type ====
    message("Detecting UMI type")
    if (any(str_detect(bcbio_nextgen_commands,
                       "work/umis/harvard-indrop-v3.json"))) {
        umi_type <- "indrop_v3"
    }


    # Sample metadata ====
    # Map the sample names to the UMIs
    if (umi_type == "indrop_v3") {
        sample_metadata <- .indrop_metadata(metadata_file, sample_dirs)
    }

    # Check that sample directories match the metadata
    if (!all(names(sample_dirs) %in% sample_metadata[["sample_barcode"]])) {
        stop("Sample name mismatch between directories and metadata")
    }


    # Sparse counts ====
    # The pipeline outputs transcript-level counts. Let's store these inside the
    # bcbioSCDataSet but outside of the SummarizedExperiment, where we will
    # instead slot the gene-level counts.
    tx_sparse_counts <- .sparse_counts(file.path(project_dir, "tagcounts.mtx"))
    sparse_counts <- .sparse_counts_tx2gene(tx_sparse_counts, genome_build)


    # Cellular barcodes ====
    cellular_barcodes <- .cellular_barcodes(sample_dirs)


    # Row data ====
    annotable <- .annotable(genome_build)


    # Column data ====
    metrics <- .metrics(sparse_counts, annotable)


    # Metadata ====
    metadata <- SimpleList(
        upload_dir = upload_dir,
        project_dir = project_dir,
        sample_dirs = sample_dirs,
        sample_metadata = sample_metadata,
        data_versions = data_versions,
        programs = programs,
        bcbio_nextgen = bcbio_nextgen,
        bcbio_nextgen_commands = bcbio_nextgen_commands,
        run_date = run_date,
        load_date = Sys.Date(),
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
