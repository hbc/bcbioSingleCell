#' Load bcbio Single-Cell RNA-Seq Run
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @param uploadDir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param sampleMetadataFile Sample barcode metadata file.
#' @param interestingGroups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param ensemblVersion Ensembl release version. Defaults to current, and does
#'   not typically need to be user-defined. This parameter can be useful for
#'   matching Ensembl annotations against an outdated bcbio annotation build.
#' @param gtfFile *Optional*. Gene transfer format (GTF) file, which will be
#'   used for transcript-to-gene (`tx2gene`) and gene-to-symbol (`gene2symbol`)
#'   annotation mappings.
#' @param wellMetadataFile *Optional*. Well identifier metadata file.
#' @param ... Additional arguments, to be stashed in the [metadata()] slot.
#'
#' @return [bcbioSingleCell].
#' @export
loadSingleCellRun <- function(
    uploadDir,
    sampleMetadataFile,
    interestingGroups = "sampleName",
    ensemblVersion = "current",
    gtfFile = NULL,
    wellMetadataFile = NULL,
    ...) {
    # Initial run setup ====
    pipeline <- "bcbio"
    if (!dir.exists(uploadDir)) {
        stop("Final upload directory does not exist", call. = FALSE)
    }
    uploadDir <- normalizePath(uploadDir)
    # Check for dated project summary directory
    detectProjectDir <-
        list.dirs(uploadDir,
                  full.names = FALSE,
                  recursive = FALSE) %>%
        str_detect("^\\d{4}-\\d{2}-\\d{2}_[^/]+$") %>%
        any()
    if (!detectProjectDir) {
        stop("Failed to locate bcbio project summary directory", call. = FALSE)
    }
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sample metadata ====
    sampleMetadataFile <- normalizePath(sampleMetadataFile)
    sampleMetadata <- .readSampleMetadataFile(
        sampleMetadataFile, sampleDirs, pipeline)

    # Check to ensure interesting groups are defined
    if (!all(interestingGroups %in% colnames(sampleMetadata))) {
        stop("Interesting groups missing in sample metadata")
    }

    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (nrow(sampleMetadata) < length(sampleDirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% rownames(sampleMetadata)]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Project directory ====
    projectDir <- dir(uploadDir,
                      pattern = projectDirPattern,
                      full.names = FALSE,
                      recursive = FALSE)
    if (length(projectDir) != 1L) {
        stop("Uncertain about project directory location")
    }
    message(projectDir)
    match <- str_match(projectDir, projectDirPattern)
    runDate <- match[[2L]] %>%
        as.Date()
    template <- match[[3L]]
    projectDir <- file.path(uploadDir, projectDir)

    # Log files ====
    message("Reading log files")
    bcbioLog <- .logFile(
        file.path(projectDir, "bcbio-nextgen.log"))
    bcbioCommandsLog <- .logFile(
        file.path(projectDir, "bcbio-nextgen-commands.log"))

    # Cellular barcode cutoff ====
    cellularBarcodeCutoffPattern <- "--cb_cutoff (\\d+)"
    cellularBarcodeCutoff <-
        str_match(bcbioCommandsLog, cellularBarcodeCutoffPattern) %>%
        .[, 2L] %>%
        na.omit() %>%
        unique() %>%
        as.numeric()

    # Data versions and programs ====
    dataVersions <- .dataVersions(projectDir)
    programs <- .programs(projectDir)
    if (!is.null(dataVersions)) {
        genomeBuild <- dataVersions %>%
            dplyr::filter(.data[["resource"]] == "transcripts") %>%
            pull("genome")
    } else {
        # Data versions aren't saved when using a custom FASTA
        # Remove this in a future update
        genomePattern <- "work/rapmap/[^/]+/quasiindex/(\\b[A-Za-z0-9]+\\b)"
        if (any(str_detect(bcbioCommandsLog, genomePattern))) {
            genomeBuild <-
                str_match(bcbioCommandsLog, genomePattern) %>%
                .[, 2L] %>%
                na.omit() %>%
                unique()
        } else {
            stop("Genome detection from bcbio commands failed")
        }
    }
    if (length(genomeBuild) > 1L) {
        stop("Multiple genomes detected -- not supported")
    }
    message(paste("Genome build:", genomeBuild))
    organism <- detectOrganism(genomeBuild)
    message(paste("Organism:", organism))

    # Molecular barcode (UMI) type ====
    umiPattern <- "/umis/([a-z0-9\\-]+)\\.json"
    if (any(str_detect(bcbioCommandsLog, umiPattern))) {
        umiType <- str_match(bcbioCommandsLog, umiPattern) %>%
            .[, 2L] %>%
            na.omit() %>%
            unique() %>%
            str_replace("-transform", "")
        message(paste("UMI type:", umiType))
    } else {
        stop("Failed to detect UMI type from JSON file")
    }

    # Well metadata ====
    if (!is.null(wellMetadataFile)) {
        wellMetadataFile <- normalizePath(wellMetadataFile)
        wellMetadata <- readFileByExtension(wellMetadataFile)
    } else {
        wellMetadata <- NULL
    }

    # tx2gene and gene2symbol annotations ====
    if (!is.null(gtfFile)) {
        gtfFile <- normalizePath(gtfFile)
        gtf <- readGTF(gtfFile)
        tx2gene <- tx2geneFromGTF(gtf)
        gene2symbol <- gene2symbolFromGTF(gtf)
    } else {
        warning(paste(
            "GTF file matching transcriptome FASTA is advised.",
            "Using tx2gene mappings from Ensembl as a fallback."
        ))
        tx2gene <- annotable(
            genomeBuild,
            format = "tx2gene",
            release = ensemblVersion)
        gene2symbol <- annotable(
            genomeBuild,
            format = "gene2symbol",
            release = ensemblVersion)
    }

    # Cellular barcodes ====
    cellularBarcodes <- .cellularBarcodesList(sampleDirs)

    # Row data =================================================================
    annotable <- annotable(genomeBuild)

    # Assays ===================================================================
    message("Reading counts")
    # Migrate this to `mapply()` method in future update
    sparseList <- pblapply(seq_along(sampleDirs), function(a) {
        txlevel <- .readSparseCounts(sampleDirs[a], pipeline = pipeline)
        # Transcript-level to gene-level counts
        genelevel <- .sparseCountsTx2Gene(txlevel, tx2gene)
        # Pre-filter using cellular barcode summary metrics
        metrics <- calculateMetrics(genelevel, annotable, prefilter = TRUE)
        genelevel[, rownames(metrics)]
    }) %>%
        setNames(names(sampleDirs))
    sparseCounts <- do.call(Matrix::cBind, sparseList)

    # Column data ==============================================================
    metrics <- calculateMetrics(sparseCounts, annotable)
    # Add reads per cellular barcode to metrics
    cellularBarcodesTibble <- .bindCellularBarcodes(cellularBarcodes) %>%
        mutate(cellularBarcode = NULL,
               sampleID = NULL)
    metrics <- metrics %>%
        as.data.frame() %>%
        rownames_to_column("cellID") %>%
        left_join(cellularBarcodesTibble, by = "cellID") %>%
        dplyr::select("nCount", everything()) %>%
        column_to_rownames("cellID")

    # Metadata =================================================================
    if (str_detect(umiType, "indrop")) {
        multiplexedFASTQ <- TRUE
    } else {
        multiplexedFASTQ <- FALSE
    }

    metadata <- list(
        version = packageVersion("bcbioSingleCell"),
        pipeline = pipeline,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = sampleMetadataFile,
        sampleMetadata = sampleMetadata,
        interestingGroups = interestingGroups,
        organism = organism,
        genomeBuild = genomeBuild,
        ensemblVersion = ensemblVersion,
        gtfFile = gtfFile,
        annotable = annotable,
        gene2symbol = gene2symbol,
        umiType = umiType,
        allSamples = allSamples,
        multiplexedFASTQ = multiplexedFASTQ,
        # bcbio pipeline-specific
        projectDir = projectDir,
        template = template,
        runDate = runDate,
        wellMetadataFile = wellMetadataFile,
        wellMetadata = wellMetadata,
        tx2gene = tx2gene,
        dataVersions = dataVersions,
        programs = programs,
        bcbioLog = bcbioLog,
        bcbioCommandsLog = bcbioCommandsLog,
        cellularBarcodeCutoff = cellularBarcodeCutoff
    )
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
        metadata <- c(metadata, dots)
    }

    # Return `bcbioSingleCell` object ==========================================
    sce <- .SingleCellExperiment(
        assays = list(assay = sparseCounts),
        colData = metrics,
        rowData = annotable,
        metadata = metadata
    )
    bcb <- new("bcbioSingleCell", sce)
    # Keep these in the bcbio slot because they contain filtered cellular
    # barcodes not present in the main assay matrix.
    bcbio(bcb, "cellularBarcodes") <- cellularBarcodes
    bcb
}
