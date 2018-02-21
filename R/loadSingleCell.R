#' Load bcbio Single-Cell RNA-Seq Data
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump annotable camel detectOrganism gene2symbolFromGTF
#'   readGTF readYAML tx2geneFromGTF
#' @importFrom bcbioBase prepareSummarizedExperiment readDataVersions
#'   readLogFile readProgramVersions readSampleMetadataFile sampleYAMLMetadata
#' @importFrom Matrix cBind
#' @importFrom pbapply pblapply
#' @importFrom stats na.omit setNames
#' @importFrom stringr str_extract str_match
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#' @param uploadDir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param sampleMetadataFile Sample barcode metadata file. Optional for runs
#'   with demultiplixed index barcodes (e.g. SureCell), but otherwise required
#'   for runs with multipliexed FASTQs containing multiple index barcodes (e.g.
#'   inDrop). Consult the GitHub repo for examples and additional information.
#' @param interestingGroups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param prefilter Prefilter counts prior to quality control analysis.
#' @param gtfFile *Optional but recommended*. GTF (Gene Transfer Format) file,
#'   which will be used for gene-to-symbol (`gene2symbol`) and
#'   transcript-to-gene (`tx2gene`) annotation mappings.
#' @param annotable User-defined gene annotations (a.k.a.
#'   "annotable"), which will be slotted into [rowData()]. Typically this should
#'   be left set to `TRUE`. By default, the function will automatically generate
#'   an annotable from the annotations available on Ensembl. If set `FALSE`,
#'   then [rowData()] inside the resulting [bcbioSingleCell] object will be left
#'   empty. This is recommended for projects dealing with genes or transcripts
#'   that are poorly annotated. Additionally, a manually constructed annotable
#'   can be passed in as a [data.frame], but this isn't generally recommended.
#' @param organism *Optional.* Organism name. Use the full latin name (e.g.
#'   "Homo sapiens"), since this will be input downstream to
#'   AnnotationHub/ensembldb. If set, this genome must be supported on Ensembl.
#'   Normally this can be left `NULL`, and the function will attempt to detect
#'   the organism automatically using [detectOrganism()].
#' @param ensemblVersion *Optional.* Ensembl release version. If `NULL`,
#'   defaults to current release, and does not typically need to be
#'   user-defined. This parameter can be useful for matching Ensembl annotations
#'   against an outdated bcbio annotation build.
#' @param genomeBuild *Optional.* Genome build. Normally this can be left `NULL`
#'   and the build will be detected from the bcbio run data. This can be set
#'   manually (e.g. "hg19" for the older *Homo sapiens* reference genome). Note
#'   that this must match the genome build identifier on Ensembl for annotations
#'   to download correctly.
#' @param ... Additional arguments, to be stashed in the [metadata()] slot.
#'
#' @return [bcbioSingleCell].
#' @export
#'
#' @examples
#' # Homo sapiens
#' extdataDir <- system.file("extdata", package = "bcbioSingleCell")
#' uploadDir <- file.path(extdataDir, "harvard_indrop_v3")
#' sampleMetadataFile <- file.path(extdataDir, "harvard_indrop_v3.xlsx")
#' bcb <- loadSingleCell(
#'     uploadDir = uploadDir,
#'     sampleMetadataFile = sampleMetadataFile)
#' print(bcb)
loadSingleCell <- function(
    uploadDir,
    sampleMetadataFile = NULL,
    interestingGroups = "sampleName",
    gtfFile = NULL,
    annotable = TRUE,
    organism = NULL,
    ensemblVersion = NULL,
    genomeBuild = NULL,
    prefilter = TRUE,
    ...) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir)
    assertIsAStringOrNULL(sampleMetadataFile)
    assert_is_character(interestingGroups)
    assertIsAStringOrNULL(gtfFile)
    if (is_a_string(gtfFile)) {
        assert_all_are_existing_files(gtfFile)
    }
    assert_is_any_of(annotable, c("data.frame", "logical", "NULL"))
    if (is.data.frame(annotable)) {
        assertIsAnnotable(annotable)
    }
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblVersion)
    assertIsAStringOrNULL(genomeBuild)
    assert_is_a_bool(prefilter)

    pipeline <- "bcbio"

    # Directory paths ==========================================================
    projectDir <- dir(
        uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE)
    assert_is_of_length(projectDir, 1L)
    inform(projectDir)
    match <- str_match(projectDir, projectDirPattern)
    runDate <- as.Date(match[[2L]])
    template <- match[[3L]]
    projectDir <- file.path(uploadDir, projectDir)
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sequencing lanes =========================================================
    if (any(grepl(lanePattern, sampleDirs))) {
        lanes <- str_match(names(sampleDirs), lanePattern) %>%
            .[, 2L] %>%
            unique() %>%
            length()
        inform(paste(
            lanes, "sequencing lane detected", "(technical replicates)"
        ))
    } else {
        lanes <- 1L
    }

    # Project summary YAML =====================================================
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- readYAML(yamlFile)

    # Log files ================================================================
    inform("Reading log files")
    bcbioLogFile <- list.files(
        path = projectDir,
        pattern = "bcbio-nextgen.log",
        full.names = TRUE,
        recursive = FALSE)
    assert_is_a_string(bcbioLogFile)
    inform(basename(bcbioLogFile))
    bcbioLog <- readLogFile(bcbioLogFile)
    assert_is_character(bcbioLog)

    bcbioCommandsLogFile <- list.files(
        path = projectDir,
        pattern = "bcbio-nextgen-commands.log",
        full.names = TRUE,
        recursive = FALSE)
    assert_is_a_string(bcbioCommandsLogFile)
    inform(basename(bcbioCommandsLogFile))
    bcbioCommandsLog <- readLogFile(bcbioCommandsLogFile)
    assert_is_character(bcbioCommandsLog)

    # Cellular barcode cutoff
    cellularBarcodeCutoffPattern <- "--cb_cutoff (\\d+)"
    assert_any_are_matching_regex(
        x = bcbioCommandsLog,
        pattern = cellularBarcodeCutoffPattern)
    match <- str_match(
        string = bcbioCommandsLog,
        pattern = cellularBarcodeCutoffPattern)
    cellularBarcodeCutoff <- match %>%
        .[, 2L] %>%
        na.omit() %>%
        unique() %>%
        as.integer()
    assert_is_an_integer(cellularBarcodeCutoff)

    inform(paste(
        cellularBarcodeCutoff,
        "reads per cellular barcode cutoff detected"
    ))

    # Detect MatrixMarket output at transcript or gene level
    genemapPattern <- "--genemap (.+)-tx2gene.tsv"
    if (any(grepl(genemapPattern, bcbioCommandsLog))) {
        countsLevel <- "gene"
    } else {
        countsLevel <- "transcript"
    }

    # Data versions and programs ===============================================
    inform("Reading data and program versions")
    dataVersions <- readDataVersions(
        file.path(projectDir, "data_versions.csv"))
    programs <- readProgramVersions(
        file.path(projectDir, "programs.txt"))

    # Detect genome build
    if (!is_a_string(genomeBuild)) {
        if (is.data.frame(dataVersions)) {
            genomeBuild <- dataVersions %>%
                .[.[["resource"]] == "transcripts", "genome", drop = TRUE]
        } else {
            # Data versions aren't saved when using a custom FASTA
            genomePattern <- "work/rapmap/[^/]+/quasiindex/(\\b[A-Za-z0-9]+\\b)"
            assert_any_are_matching_regex(bcbioCommandsLog, genomePattern)
            genomeBuild <- str_match(
                bcbioCommandsLog, genomePattern) %>%
                .[, 2L] %>%
                na.omit() %>%
                unique()
        }
    }
    assert_is_a_string(genomeBuild)

    # Detect organism
    if (!is_a_string(organism)) {
        assert_is_a_string(genomeBuild)
        organism <- detectOrganism(genomeBuild)
    }
    assert_is_a_string(organism)

    inform(paste(
        paste("Organism:", organism),
        paste("Genome build:", genomeBuild),
        sep = "\n"
    ))

    # Molecular barcode (UMI) type =============================================
    umiPattern <- "/umis/([a-z0-9\\-]+)\\.json"
    assert_any_are_matching_regex(bcbioCommandsLog, umiPattern)
    umiType <- str_match(bcbioCommandsLog, umiPattern) %>%
        .[, 2L] %>%
        na.omit() %>%
        unique() %>%
        gsub(
            x = .,
            pattern = "-transform",
            replacement = "")
    assert_is_a_string(umiType)
    inform(paste("UMI type:", umiType))

    # Sample metadata ==========================================================
    if (is_a_string(sampleMetadataFile)) {
        sampleMetadata <- readSampleMetadataFile(sampleMetadataFile)
    } else {
        assert_all_are_matching_regex(umiType, "indrop")
        sampleMetadata <- sampleYAMLMetadata(yaml)
    }

    if ("sequence" %in% colnames(sampleMetadata)) {
        sampleDirSequence <- str_extract(names(sampleDirs), "[ACGT]+$")
        if (identical(
            sort(sampleDirSequence),
            sort(as.character(sampleMetadata[["sequence"]]))
        )) {
            warn(paste(
                "It appears that the reverse complement sequence of the",
                "i5 index barcode(s) was input into the sample metadata",
                "`sequence` column. bcbio outputs the revcomp into the",
                "sample directories, but the forward sequence should be",
                "used in the R package."
            ))
        }
    }
    assert_are_identical(rownames(sampleMetadata), names(sampleDirs))

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    assertFormalInterestingGroups(sampleMetadata, interestingGroups)

    # Subset sample directories by metadata ====================================
    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (nrow(sampleMetadata) < length(sampleDirs)) {
        inform("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% rownames(sampleMetadata)]
        inform(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Gene annotations =========================================================
    # Ensembl annotations (gene annotable)
    if (isTRUE(annotable)) {
        annotable <- annotable(
            organism,
            genomeBuild = genomeBuild,
            release = ensemblVersion,
            uniqueSymbol = FALSE)
    } else if (is.data.frame(annotable)) {
        annotable <- annotable(annotable)
    } else {
        warn("Loading run without gene annotations")
        annotable <- NULL
    }

    # GTF annotations
    if (is_a_string(gtfFile)) {
        gtf <- readGTF(gtfFile)
    } else {
        gtf <- NULL
    }

    # Transcript-to-gene mappings
    if (countsLevel == "transcript") {
        assert_is_data.frame(gtf)
        tx2gene <- tx2geneFromGTF(gtf)
    } else {
        # Not applicable to gene-level bcbio output
        tx2gene <- NA
    }

    # Gene-to-symbol mappings
    if (is_a_string(gtfFile)) {
        gene2symbol <- gene2symbolFromGTF(gtf)
    } else if (is.data.frame(annotable)) {
        gene2symbol <- annotable[, c("ensgene", "symbol")]
    } else {
        abort("Loading run without gene-to-symbol mappings (not recommended)")
        gene2symbol <- NULL
    }

    # Cellular barcode distributions ===========================================
    cbList <- .cellularBarcodesList(sampleDirs)
    cbData <- .bindCellularBarcodes(cbList)

    # Counts ===================================================================
    inform(paste("Reading counts at", countsLevel, "level"))
    sparseCountsList <- .sparseCountsList(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        umiType = umiType)
    # Combine the individual per-sample transcript-level sparse matrices into a
    # single sparse matrix
    counts <- do.call(Matrix::cBind, sparseCountsList)
    # Convert counts from transcript-level to gene-level, if necessary
    if (countsLevel == "transcript") {
        counts <- .transcriptToGeneLevelCounts(counts, tx2gene)
    }

    # Metrics ==================================================================
    metrics <- calculateMetrics(
        object = counts,
        annotable = annotable,
        prefilter = prefilter)
    # Bind the `nCount` column to the metrics
    cbPass <- cbData[rownames(metrics), "nCount", drop = FALSE]
    metrics <- cbind(metrics, cbPass)

    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the cells that passed prefiltering
        counts <- counts[, rownames(metrics), drop = FALSE]
    }

    # Cell to sample mappings ==================================================
    cell2sample <- mapCellsToSamples(
        cells = rownames(metrics),
        samples = rownames(sampleMetadata))

    # Metadata =================================================================
    metadata <- list(
        version = packageVersion,
        pipeline = pipeline,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = sampleMetadataFile,
        sampleMetadata = sampleMetadata,
        interestingGroups = interestingGroups,
        cell2sample = cell2sample,
        organism = organism,
        genomeBuild = genomeBuild,
        ensemblVersion = ensemblVersion,
        gtfFile = gtfFile,
        annotable = annotable,
        gene2symbol = gene2symbol,
        umiType = umiType,
        allSamples = allSamples,
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
        cellularBarcodeCutoff = cellularBarcodeCutoff)
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
        metadata <- c(metadata, dots)
    }

    # Return ===================================================================
    # Use an internal `SummarizedExperiment()` function call to handle rowname
    # mismatches with the annotable. This can happen when newer Ensembl
    # annotations are requested than those used for count alignment, or when
    # we pass in FASTA spike-ins (e.g. EGFP).
    se <- prepareSummarizedExperiment(
        assays = list(assay = counts),
        rowData = annotable,
        colData = metrics,
        metadata = metadata)
    bcbio <- list(cellularBarcodes = cbList)
    new("bcbioSingleCell", se, bcbio = as(bcbio, "SimpleList"))
}
