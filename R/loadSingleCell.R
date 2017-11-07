#' Load bcbio Single-Cell RNA-Seq Data
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump annotable camel detectOrganism gene2symbolFromGTF
#'   prepareSummarizedExperiment readDataVersions readGTF readLogFile
#'   readProgramVersions readSampleMetadataFile readYAML sampleYAMLMetadata
#'   tx2geneFromGTF
#' @importFrom dplyr everything filter left_join mutate pull select
#' @importFrom Matrix cBind
#' @importFrom pbapply pblapply
#' @importFrom stats na.omit setNames
#' @importFrom stringr str_match
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom utils packageVersion
#'
#' @param uploadDir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param interestingGroups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param sampleMetadataFile Sample barcode metadata file.
#' @param gtfFile *Optional*. GTF (Gene Transfer Format) file, which will be
#'   used for transcript-to-gene (`tx2gene`) and gene-to-symbol (`gene2symbol`)
#'   annotation mappings.
#' @param prefilter Prefilter counts prior to quality control analysis.
#' @param annotable *Optional*. User-defined gene annotations (a.k.a.
#'   "annotable"), which will be slotted into [rowData()]. Typically this should
#'   be left undefined. By default, the function will automatically generate an
#'   annotable from the annotations available on Ensembl. If set `NULL`, then
#'   [rowData()] inside the resulting [bcbioSingleCell] object will be left
#'   empty. This is recommended for projects dealing with genes or transcripts
#'   that are poorly annotated.
#' @param ensemblVersion *Optional*. Ensembl release version. If `NULL`,
#'   defaults to current release, and does not typically need to be
#'   user-defined. This parameter can be useful for matching Ensembl annotations
#'   against an outdated bcbio annotation build.
#' @param ... Additional arguments, to be stashed in the [metadata()] slot.
#'
#' @return [bcbioSingleCell].
#' @export
#'
#' @examples
#' \dontrun{
#' # Mus musculus
#' # Run with Ensembl 88 transcriptome FASTA and GTF files
#' bcb <- loadSingleCell(
#'     uploadDir = file.path("indrop_rnaseq", "final"),
#'     interestingGroups = c("genotype", "treatment"),
#'     sampleMetadataFile = file.path("meta", "sample_metadata.xlsx"),
#'     gtfFile = file.path(
#'         "annotations",
#'         "Mus_musculus.GRCm38.88.chr_patch_hapl_scaff.gtf.gz"),
#'     ensemblVersion = 88)
#' }
loadSingleCell <- function(
    uploadDir,
    interestingGroups = "sampleName",
    sampleMetadataFile = NULL,
    gtfFile = NULL,
    prefilter = TRUE,
    annotable,
    ensemblVersion = NULL,
    ...) {
    pipeline <- "bcbio"

    # Directory paths ==========================================================
    # Check connection to final upload directory
    if (!dir.exists(uploadDir)) {
        stop("Final upload directory does not exist", call. = FALSE)
    }
    uploadDir <- normalizePath(uploadDir)
    projectDir <- dir(
        uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE)
    if (length(projectDir) != 1) {
        stop("Uncertain about project directory location", call. = FALSE)
    }
    message(projectDir)
    match <- str_match(projectDir, projectDirPattern)
    runDate <- as.Date(match[[2]])
    template <- match[[3]]
    projectDir <- file.path(uploadDir, projectDir)
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sequencing lanes =========================================================
    if (any(grepl(x = sampleDirs, pattern = lanePattern))) {
        lanes <- str_match(names(sampleDirs), lanePattern) %>%
            .[, 2] %>%
            unique() %>%
            length()
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"))
    } else {
        lanes <- 1
    }

    # Project summary YAML =====================================================
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- readYAML(yamlFile)

    # Log files ================================================================
    message("Reading log files")
    bcbioLog <- readLogFile(
        file.path(projectDir, "bcbio-nextgen.log"))
    bcbioCommandsLog <- readLogFile(
        file.path(projectDir, "bcbio-nextgen-commands.log"))

    # Cellular barcode cutoff
    cellularBarcodeCutoffPattern <- "--cb_cutoff (\\d+)"
    cellularBarcodeCutoff <-
        str_match(bcbioCommandsLog, cellularBarcodeCutoffPattern) %>%
        .[, 2] %>%
        na.omit() %>%
        unique() %>%
        as.numeric()

    # Detect MatrixMarket output at transcript or gene level
    # This grep pattern may not be strict enough against the file path
    genemapPattern <- "--genemap (.+)-tx2gene.tsv"
    if (any(grepl(x = bcbioCommandsLog, pattern = genemapPattern))) {
        countsLevel <- "gene"
    } else {
        countsLevel <- "transcript"
    }

    # Data versions and programs ===============================================
    dataVersions <- readDataVersions(
        file.path(projectDir, "data_versions.csv"))
    programs <- readProgramVersions(
        file.path(projectDir, "programs.txt"))
    if (!is.null(dataVersions)) {
        genomeBuild <- dataVersions %>%
            filter(.data[["resource"]] == "transcripts") %>%
            pull("genome")
    } else {
        # Data versions aren't saved when using a custom FASTA
        # Remove this in a future update
        genomePattern <- "work/rapmap/[^/]+/quasiindex/(\\b[A-Za-z0-9]+\\b)"
        if (any(grepl(x = bcbioCommandsLog, pattern = genomePattern))) {
            genomeBuild <-
                str_match(bcbioCommandsLog, genomePattern) %>%
                .[, 2] %>%
                na.omit() %>%
                unique()
        } else {
            stop("Genome detection from bcbio commands failed", call. = FALSE)
        }
    }
    if (length(genomeBuild) > 1) {
        stop("Multiple genomes detected", call. = FALSE)
    }
    organism <- detectOrganism(genomeBuild)
    message(paste0("Genome: ", organism, " (", genomeBuild, ")"))

    # Molecular barcode (UMI) type =============================================
    umiPattern <- "/umis/([a-z0-9\\-]+)\\.json"
    if (any(grepl(x = bcbioCommandsLog, pattern = umiPattern))) {
        umiType <- str_match(bcbioCommandsLog, umiPattern) %>%
            .[, 2] %>%
            na.omit() %>%
            unique() %>%
            gsub(x = .,
                 pattern = "-transform",
                 replacement = "")
        message(paste("UMI type:", umiType))
    } else {
        stop("Failed to detect UMI type from JSON file", call. = FALSE)
    }

    # Multiplexed FASTQ ========================================================
    # This value determines how we assign sampleIDs and downstream plot
    # appearance in the quality control analysis
    if (grepl(x = umiType, pattern = "indrop")) {
        multiplexedFASTQ <- TRUE
    } else {
        multiplexedFASTQ <- FALSE
    }

    # Sample metadata ==========================================================
    if (!is.null(sampleMetadataFile)) {
        sampleMetadataFile <- normalizePath(sampleMetadataFile)
        sampleMetadata <- readSampleMetadataFile(sampleMetadataFile)
    } else {
        if (grepl(x = umiType, pattern = "indrop")) {
            # Enforce `sampleMetadataFile` for multiplexed samples
            stop("'sampleMetadataFile' is required for multiplexed samples",
                 call. = FALSE)
        }
        sampleMetadata <- sampleYAMLMetadata(yaml)
    }
    # Check that `sampleID` matches `sampleDirs`
    if (!all(sampleMetadata[["sampleID"]] %in% names(sampleDirs))) {
        stop("Sample directory names don't match the sample metadata file",
             call. = FALSE)
    }

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    # Check to ensure interesting groups are defined
    if (!all(interestingGroups %in% colnames(sampleMetadata))) {
        stop("Interesting groups missing in sample metadata", call. = FALSE)
    }

    # Subset sample directories by metadata ====================================
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

    # Cellular barcodes ========================================================
    cellularBarcodes <- .cellularBarcodesList(sampleDirs)

    # Gene annotations =========================================================
    if (missing(annotable)) {
        annotable <- basejump::annotable(
            organism,
            genomeBuild = genomeBuild,
            release = ensemblVersion)
    } else if (is.data.frame(annotable)) {
        annotable <- annotable(annotable)
    }
    if (!is.null(gtfFile)) {
        gtfFile <- normalizePath(gtfFile)
        gtf <- readGTF(gtfFile)
        if (countsLevel == "transcript") {
            tx2gene <- tx2geneFromGTF(gtf)
        }
        gene2symbol <- gene2symbolFromGTF(gtf)
    } else {
        warning(paste(
            "GFF/GTF file matching transcriptome FASTA is advised.",
            "Using tx2gene mappings from Ensembl as a fallback."
        ), call. = FALSE)
        if (countsLevel == "transcript") {
            tx2gene <- annotable(
                organism,
                format = "tx2gene",
                release = ensemblVersion)
        }
        gene2symbol <- annotable(
            organism,
            format = "gene2symbol",
            release = ensemblVersion)
    }

    # Counts ===================================================================
    message("Reading counts")
    # Migrate this to `mapply()` method in future update
    sparseList <- pblapply(seq_along(sampleDirs), function(a) {
        .readSparseCounts(
            sampleDirs[a],
            pipeline = pipeline,
            umiType = umiType)
    }) %>%
        setNames(names(sampleDirs))
    # Combine the individual per-sample transcript-level sparse matrices into a
    # single sparse matrix
    counts <- do.call(Matrix::cBind, sparseList)
    # Convert counts from transcript-level to gene-level, if necessary
    if (countsLevel == "transcript") {
        counts <- .sparseCountsTx2Gene(counts, tx2gene)
    }

    # Metrics ==================================================================
    metrics <- calculateMetrics(
        counts,
        annotable = annotable,
        prefilter = prefilter)
    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the cells that passed prefiltering
        counts <- counts[, rownames(metrics)]
    }

    # Slot the cellular barcode metrics into colData
    # Add reads per cellular barcode
    cellularBarcodesTibble <- .bindCellularBarcodes(cellularBarcodes) %>%
        mutate(cellularBarcode = NULL,
               sampleID = NULL)
    metrics <- metrics %>%
        as.data.frame() %>%
        rownames_to_column("cellID") %>%
        left_join(cellularBarcodesTibble, by = "cellID") %>%
        select("nCount", everything()) %>%
        column_to_rownames("cellID")

    # Metadata =================================================================
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
        prefilter = prefilter,
        # bcbio pipeline-specific
        projectDir = projectDir,
        template = template,
        runDate = runDate,
        yamlFile = yamlFile,
        yaml = yaml,
        tx2gene = tx2gene,
        dataVersions = dataVersions,
        programs = programs,
        bcbioLog = bcbioLog,
        bcbioCommandsLog = bcbioCommandsLog,
        cellularBarcodeCutoff = cellularBarcodeCutoff
    )
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0) {
        metadata <- c(metadata, dots)
    }

    # SummarizedExperiment =====================================================
    # Use an internal `SummarizedExperiment()` function call to handle rowname
    # mismatches with the annotable. This can happen when newer Ensembl
    # annotations are requested than those used for count alignment, or when
    # we pass in FASTA spike-ins (e.g. EGFP).
    se <- prepareSummarizedExperiment(
        assays = list(assay = counts),
        rowData = annotable,
        colData = metrics,
        metadata = metadata
    )
    bcb <- new("bcbioSingleCell", se)
    # Keep these in the bcbio slot because they contain filtered cellular
    # barcodes not present in the main assay matrix.
    bcbio(bcb, "cellularBarcodes") <- cellularBarcodes
    bcb
}
