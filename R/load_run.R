#' Load `bcbio` Run
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
#' @param tx2gene tx2gene [data.frame], required when using a custom
#'   transcriptome FASTA file. This argument will be removed in a future update.
#' @param ... Additional arguments, passed as metadata.
#'
#' @return [bcbioSCDataSet].
#' @export
load_run <- function(
    upload_dir = "final",
    sample_metadata_file,
    well_metadata_file = NULL,
    interesting_groups = "sample_name",
    tx2gene = NULL,
    ...) {
    if (!dir.exists(upload_dir)) {
        stop("Upload directory missing")
    }


    # Initial run setup ====
    upload_dir <- normalizePath(upload_dir)
    if (!dir.exists(upload_dir)) {
        stop("Final upload directory does not exist")
    }
    pipeline <- .detect_pipeline(upload_dir)
    sample_dirs <- .sample_dirs(upload_dir, pipeline = pipeline)


    # Sample metadata ====
    sample_metadata_file <- normalizePath(sample_metadata_file)
    sample_metadata <- .sample_metadata_file(
        sample_metadata_file, sample_dirs, pipeline)

    # Check to ensure interesting groups are defined
    if (!all(interesting_groups %in% colnames(sample_metadata))) {
        stop("Interesting groups missing in sample metadata")
    }

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
        # Project directory ----
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

        # Log files ----
        message("Reading log files")
        bcbio_log <- .log_file(
            file.path(project_dir, "bcbio-nextgen.log"))
        bcbio_commands_log <- .log_file(
            file.path(project_dir, "bcbio-nextgen-commands.log"))

        # Cellular barcode cutoff ----
        cb_cutoff_pattern <- "--cb_cutoff (\\d+)"
        cb_cutoff <- str_match(bcbio_commands_log, cb_cutoff_pattern) %>%
            .[, 2L] %>%
            na.omit %>%
            unique %>%
            as.numeric

        # Data versions and programs ----
        data_versions <- .data_versions(project_dir)
        programs <- .programs(project_dir)
        if (!is.null(data_versions)) {
            genome_build <- data_versions %>%
                filter(.data[["resource"]] == "transcripts") %>%
                pull("genome")
        } else {
            # Data versions aren't saved when using a custom FASTA
            # Remove this in a future update
            genome_pattern <- "work/rapmap/[^/]+/quasiindex/([^/]+)/"
            if (any(str_detect(bcbio_commands_log, genome_pattern))) {
                genome_build <- str_match(bcbio_commands_log,
                                          genome_pattern) %>%
                    .[, 2L] %>%
                    na.omit %>%
                    unique
            } else {
                stop("Genome detection from bcbio commands failed")
            }
        }
        if (length(genome_build) > 1L) {
            stop("Multiple genomes detected -- not supported")
        }

        # Molecular barcode (UMI) type ----
        umi_pattern <- "/umis/([a-z0-9\\-]+)\\.json"
        if (any(str_detect(bcbio_commands_log, umi_pattern))) {
            umi_type <- str_match(bcbio_commands_log,
                                  umi_pattern) %>%
                .[, 2L] %>%
                na.omit %>%
                unique %>%
                str_replace("-transform", "")
        } else {
            stop("Failed to detect UMI type from JSON file")
        }

        # Well metadata ----
        if (!is.null(well_metadata_file)) {
            well_metadata_file <- normalizePath(well_metadata_file)
        }
        well_metadata <- readFileByExtension(well_metadata_file)

        # tx2gene ----
        if (is.null(tx2gene)) {
            tx2gene <- tx2gene(genome_build)
        }

        # Cellular barcodes ----
        cellular_barcodes <- .cellular_barcodes(sample_dirs)
    } else if (pipeline == "cellranger") {
        # Get genome build from sample_dirs
        genome_build <- basename(sample_dirs) %>% unique
        if (length(genome_build) > 1L) {
            stop("Multiple genomes detected -- not supported")
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
                .sparse_counts_tx2gene(sparse_counts, tx2gene)
        }
        # Pre-filter using cellular barcode summary metrics
        metrics <- .calculate_metrics(sparse_counts, annotable)
        sparse_counts[, rownames(metrics)]
    }) %>%
        set_names(names(sample_dirs))
    sparse_counts <- do.call(cBind, sparse_list)


    # Column data ====
    metrics <- .calculate_metrics(sparse_counts, annotable)
    if (pipeline == "bcbio") {
        # Add reads per cellular barcode to metrics
        cb_tbl <- .bind_cellular_barcodes(cellular_barcodes) %>%
            mutate(cellular_barcode = NULL,
                   sample_id = NULL)
        metrics <- metrics %>%
            as.data.frame %>%
            rownames_to_column %>%
            left_join(cb_tbl, by = "rowname") %>%
            tidy_select("reads", everything()) %>%
            column_to_rownames %>%
            as.matrix
    }


    # Metadata ====
    if (umi_type == "indrop") {
        multiplexed_fastq <- TRUE
    } else {
        multiplexed_fastq <- FALSE
    }
    metadata <- SimpleList(
        pipeline = pipeline,
        upload_dir = upload_dir,
        sample_dirs = sample_dirs,
        sample_metadata_file = sample_metadata_file,
        sample_metadata = sample_metadata,
        interesting_groups = interesting_groups,
        genome_build = genome_build,
        annotable = annotable,
        ensembl_version = ensemblVersion(),
        umi_type = umi_type,
        all_samples = all_samples,
        multiplexed_fastq = multiplexed_fastq)
    if (pipeline == "bcbio") {
        bcbio_metadata <- SimpleList(
            project_dir = project_dir,
            well_metadata_file = well_metadata_file,
            well_metadata = well_metadata,
            template = template,
            run_date = run_date,
            tx2gene = tx2gene,
            data_versions = data_versions,
            programs = programs,
            bcbio_log = bcbio_log,
            bcbio_commands_log = bcbio_commands_log,
            cb_cutoff = cb_cutoff)
        metadata <- c(metadata, bcbio_metadata)
    }
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
        metadata <- c(metadata, dots)
    }


    # Return ====
    se <- packageSE(
        assays = SimpleList(
            sparse_counts = sparse_counts),
        colData = metrics,
        rowData = annotable,
        metadata = metadata)
    bcb <- new("bcbioSCDataSet", se)
    if (pipeline == "bcbio") {
        bcbio(bcb, "cellular_barcodes") <- cellular_barcodes
    }
    bcb
}
