#' Load bcbio-nextgen run
#'
#' @details
#' - **cellranger**: Read [10x Genomics
#'   Chromium](https://www.10xgenomics.com/software/) cell counts from
#'   `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param sample_metadata_file **Required**. Sample barcode metadata file.
#' @param well_metadata_file *Optional*. Well identifier metadata file.
#' @param interesting_groups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param ... Additional arguments, passed as metadata.
#'
#' @return [bcbioSCDataSet].
#' @export
load_run <- function(
    upload_dir = "final",
    sample_metadata_file,
    well_metadata_file = NULL,
    interesting_groups = "sample_name",
    ...) {
    if (!dir.exists(upload_dir)) {
        stop("Upload directory missing")
    }


    # Initial run setup ====
    upload_dir <- normalizePath(upload_dir)
    pipeline <- .detect_pipeline(upload_dir)
    sample_dirs <- .sample_dirs(upload_dir, pipeline = pipeline)


    # Sample metadata ====
    sample_metadata_file <- normalizePath(sample_metadata_file)
    sample_metadata <- .sample_metadata_file(
        sample_metadata_file, sample_dirs, pipeline)

    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (length(sample_metadata[["sample_id"]]) < length(sample_dirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        all_samples <- FALSE
        sample_dirs <- sample_dirs %>%
            .[names(sample_dirs) %in% sample_metadata[["sample_id"]]]
        message(paste(length(sample_dirs), "samples matched by metadata"))
    } else {
        all_samples <- TRUE
    }


    # Pipeline-specific support prior to count loading ====
    if (pipeline == "bcbio") {
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

        # Data versions and programs ----
        data_versions <- .data_versions(project_dir)
        programs <- .programs(project_dir)
        genome_build <- data_versions %>%
            filter(.data[["resource"]] == "transcripts") %>%
            pull("genome")

        # Log files ----
        message("Reading log files")
        bcbio_nextgen <- read_lines(
            file.path(project_dir, "bcbio-nextgen.log"))
        bcbio_nextgen_commands <- read_lines(
            file.path(project_dir, "bcbio-nextgen-commands.log"))

        # Molecular barcode (UMI) type ----
        if (any(str_detect(bcbio_nextgen_commands,
                           "work/umis/harvard-indrop-v3.json"))) {
            umi_type <- "harvard-indrop-v3"
        }

        # Cellular barcodes ----
        cellular_barcodes <- .cellular_barcodes(sample_dirs)

        # Well metadata ----
        if (!is.null(well_metadata_file)) {
            well_metadata_file <- normalizePath(well_metadata_file)
        }
        well_metadata <- read_file_by_extension(well_metadata_file)
    } else if (pipeline == "cellranger") {
        # Get genome build from sample_dirs
        genome_build <- basename(sample_dirs) %>% unique
        if (length(genome_build) > 1L) {
            stop("Multiple genomes detected in cellranger samples")
        }
        umi_type <- "chromium"
    }
    message(paste("UMI type:", umi_type))


    # Row data ====
    annotable <- annotable(genome_build)


    # Read counts into sparse matrix ====
    message("Reading counts")
    sparse_list <- pblapply(seq_along(sample_dirs), function(a) {
        sparse_counts <- .sparse_counts(sample_dirs[a], pipeline = pipeline)
        if (pipeline == "bcbio") {
            # Convert transcript-level to gene-level
            sparse_counts <-
                .sparse_counts_tx2gene(sparse_counts, genome_build)
        }
        # Pre-filter by calculating metrics
        metrics <- .calculate_metrics(sparse_counts, annotable)
        sparse_counts[, rownames(metrics)]
    }) %>%
        set_names(names(sample_dirs))
    sparse_counts <- do.call(cBind, sparse_list)
    rm(sparse_list)


    # Column data ====
    metrics <- .calculate_metrics(sparse_counts, annotable)
    if (pipeline == "bcbio") {
        # Add reads per cellular barcode to bcbio-nextgen metrics
        cb_df <- .bind_cellular_barcodes(cellular_barcodes) %>%
            filter(.data[["cellular_barcode"]] %in% rownames(metrics))
        metrics <- metrics %>%
            as("data.frame") %>%
            rownames_to_column %>%
            left_join(cb_df,
                      by = c("rowname" = "cellular_barcode")) %>%
            tidy_select("reads", everything()) %>%
            column_to_rownames %>%
            as("matrix")
    }


    # Metadata ====
    metadata <- SimpleList(
        pipeline = pipeline,
        upload_dir = upload_dir,
        sample_dirs = sample_dirs,
        sample_metadata_file = sample_metadata_file,
        sample_metadata = sample_metadata,
        interesting_groups = interesting_groups,
        genome_build = genome_build,
        annotable = annotable,
        ensembl_version = annotables::ensembl_version,
        umi_type = umi_type,
        all_samples = all_samples)
    if (pipeline == "bcbio") {
        bcbio_metadata <- SimpleList(
            project_dir = project_dir,
            well_metadata_file = well_metadata_file,
            well_metadata = well_metadata,
            template = template,
            run_date = run_date,
            tx2gene = annotable(genome_build, format = "tx2gene"),
            data_versions = data_versions,
            programs = programs,
            bcbio_nextgen = bcbio_nextgen,
            bcbio_nextgen_commands = bcbio_nextgen_commands)
        metadata <- c(metadata, bcbio_metadata)
    }
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
        metadata <- c(metadata, dots)
    }


    # Return ====
    se <- .summarized_experiment(
        assays = SimpleList(
            sparse_counts = sparse_counts),
        col_data = metrics,
        row_data = annotable,
        metadata = metadata)
    bcb <- new("bcbioSCDataSet", se)
    if (pipeline == "bcbio") {
        bcbio(bcb, "cellular_barcodes") <- cellular_barcodes
    }
    bcb
}
