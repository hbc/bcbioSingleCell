#' Load bcbio Single-Cell RNA-Seq Data
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump camel detectOrganism ensembl readYAML
#' @importFrom bcbioBase prepareSummarizedExperiment readDataVersions
#'   readLogFile readProgramVersions readSampleMetadataFile sampleYAMLMetadata
#' @importFrom GenomicFeatures exonsBy makeTxDbFromGFF
#' @importFrom Matrix cBind
#' @importFrom pbapply pblapply
#' @importFrom SingleCellExperiment isSpike<- SingleCellExperiment
#' @importFrom stats setNames
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
#' @param organism *Optional.* Organism name. Use the full latin name (e.g.
#'   "Homo sapiens"), since this will be input downstream to
#'   AnnotationHub/ensembldb. If set, this genome must be supported on Ensembl.
#'   Normally this can be left `NULL`, and the function will attempt to detect
#'   the organism automatically using [detectOrganism()].
#' @param ensemblRelease *Optional.* Ensembl release version. If `NULL`,
#'   defaults to current release, and does not typically need to be
#'   user-defined. This parameter can be useful for matching Ensembl annotations
#'   against an outdated bcbio annotation build.
#' @param genomeBuild *Optional.* Genome build. Normally this can be left `NULL`
#'   and the build will be detected from the bcbio run data. This can be set
#'   manually (e.g. "hg19" for the older *Homo sapiens* reference genome). Note
#'   that this must match the genome build identifier on Ensembl for annotations
#'   to download correctly.
#' @param isSpike *Optional.* Gene names corresponding to FASTA spike-in
#'   sequences (e.g. ERCCs, EGFP, TDTOMATO).
#' @param gffFile *Optional.* By default, we recommend leaving this `NULL` for
#'   genomes that are supported on Ensembl. In this case, the row annotations
#'   ([rowRanges()]) will be obtained automatically from Ensembl using
#'   AnnotationHub and ensembldb. Internally, they are stored as genomic ranges
#'   (`GRanges`). Additionally, these annotations are accessible as a
#'   `data.frame` via the [rowData()] function. For a genome that is not
#'   supported on Ensembl or AnnotationHub, a GFF/GTF (General Feature Format)
#'   file is required. Generally, we recommend using a GTF (GFFv2) file here
#'   over a GFF3 file, although both formats are supported. The function will
#'   internally generate a `TxDb` containing transcript-to-gene mappings and
#'   construct a `GRanges` object containing the genomic ranges ([rowRanges()]).
#' @param prefilter Prefilter counts prior to quality control analysis.
#' @param ... Additional arguments, to be stashed in the [metadata()] slot.
#'
#' @return `bcbioSingleCell`.
#' @export
#'
#' @examples
#' extdataDir <- system.file("extdata", package = "bcbioSingleCell")
#' uploadDir <- file.path(extdataDir, "harvard_indrop_v3")
#' sampleMetadataFile <- file.path(extdataDir, "harvard_indrop_v3.xlsx")
#' loadSingleCell(
#'     uploadDir = uploadDir,
#'     sampleMetadataFile = sampleMetadataFile
#' )
loadSingleCell <- function(
    uploadDir,
    sampleMetadataFile = NULL,
    interestingGroups = "sampleName",
    organism = NULL,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    isSpike = NULL,
    gffFile = NULL,
    prefilter = TRUE,
    ...
) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir)
    assertIsAStringOrNULL(sampleMetadataFile)
    assert_is_character(interestingGroups)
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsAStringOrNULL(genomeBuild)
    assertIsCharacterOrNULL(isSpike)
    assertIsAStringOrNULL(gffFile)
    if (is_a_string(gffFile)) {
        assert_all_are_existing_files(gffFile)
    }
    assert_is_a_bool(prefilter)
    dots <- list(...)
    pipeline <- "bcbio"

    # Legacy arguments =========================================================
    call <- match.call(expand.dots = TRUE)
    # annotable
    if ("annotable" %in% names(call)) {
        abort("Use `gffFile` instead of `annotable`")
    }
    # ensemblVersion
    if ("ensemblVersion" %in% names(call)) {
        warn("Use `ensemblRelease` instead of `ensemblVersion`")
        ensemblRelease <- call[["ensemblVersion"]]
        dots[["ensemblVersion"]] <- NULL
    }
    # gtfFile
    if ("gtfFile" %in% names(call)) {
        warn("Use `gffFile` instead of `gtfFile`")
        gffFile <- call[["gtfFile"]]
        dots[["gtfFile"]] <- NULL
    }
    dots <- Filter(Negate(is.null), dots)

    # Directory paths ==========================================================
    projectDir <- dir(
        uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE
    )
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
        recursive = FALSE
    )
    assert_is_a_string(bcbioLogFile)
    inform(basename(bcbioLogFile))
    bcbioLog <- readLogFile(bcbioLogFile)
    assert_is_character(bcbioLog)

    bcbioCommandsLogFile <- list.files(
        path = projectDir,
        pattern = "bcbio-nextgen-commands.log",
        full.names = TRUE,
        recursive = FALSE
    )
    assert_is_a_string(bcbioCommandsLogFile)
    inform(basename(bcbioCommandsLogFile))
    bcbioCommandsLog <- readLogFile(bcbioCommandsLogFile)
    assert_is_character(bcbioCommandsLog)

    # Cellular barcode cutoff
    cellularBarcodeCutoffPattern <- "--cb_cutoff (\\d+)"
    assert_any_are_matching_regex(
        x = bcbioCommandsLog,
        pattern = cellularBarcodeCutoffPattern
    )
    match <- str_match(
        string = bcbioCommandsLog,
        pattern = cellularBarcodeCutoffPattern
    )
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
        level <- "gene"
    } else {
        level <- "transcript"
    }

    # Data and program versions ================================================
    inform("Reading data and program versions")
    dataVersions <- readDataVersions(
        file.path(projectDir, "data_versions.csv")
    )
    programVersions <- readProgramVersions(
        file.path(projectDir, "programs.txt")
    )

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
                string = bcbioCommandsLog,
                pattern = genomePattern
            ) %>%
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
            pattern = "-transform",
            replacement = "",
            x = .
        )
    assert_is_a_string(umiType)
    inform(paste("UMI type:", umiType))

    # Sample metadata ==========================================================
    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleMetadataFile(sampleMetadataFile)
    } else {
        assert_all_are_matching_regex(umiType, "indrop")
        sampleData <- sampleYAMLMetadata(yaml)
    }
    if ("sequence" %in% colnames(sampleData)) {
        sampleDirSequence <- str_extract(names(sampleDirs), "[ACGT]+$")
        if (identical(
            sort(sampleDirSequence),
            sort(as.character(sampleData[["sequence"]]))
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
    assert_are_identical(rownames(sampleData), names(sampleDirs))

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    assertFormalInterestingGroups(sampleData, interestingGroups)

    # Subset sample directories by metadata ====================================
    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (nrow(sampleData) < length(sampleDirs)) {
        inform("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% rownames(sampleData)]
        inform(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Unfiltered cellular barcode distributions ================================
    cbList <- .cellularBarcodesList(sampleDirs)
    cbData <- .bindCellularBarcodes(cbList)

    # Assays ===================================================================
    inform(paste("Reading counts at", level, "level"))
    sparseCountsList <- .sparseCountsList(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        umiType = umiType
    )
    counts <- do.call(Matrix::cBind, sparseCountsList)

    # Row data =================================================================
    ahMeta <- NULL
    txdb <- NULL
    if (is_a_string(gffFile)) {
        txdb <- makeTxDbFromGFF(gffFile)
        # rowRanges <- exonsBy(txdb, by = "gene")
        rowRanges <- genes(txdb)
        # Transcript-to-gene mappings
        if (level == "transcripts") {
            transcripts <- transcripts(txdb, columns = c("tx_name", "gene_id"))
            tx2gene <- mcols(transcripts) %>%
                as.data.frame() %>%
                set_colnames(c("txID", "geneID")) %>%
                set_rownames(.[["txID"]])
        }
    } else {
        # ah = AnnotationHub
        ah <- ensembl(
            organism = organism,
            format = "genes",
            genomeBuild = genomeBuild,
            release = ensemblRelease,
            return = "GRanges",
            metadata = TRUE
        )
        assert_is_list(ah)
        assert_are_identical(names(ah), c("data", "metadata"))
        rowRanges <- ah[["data"]]
        assert_is_all_of(rowRanges, "GRanges")
        ahMeta <- ah[["metadata"]]
        assert_all_are_matching_regex(ahMeta[["id"]], "^AH\\d+$")
        # Transcript-to-gene mappings
        if (level == "transcripts") {
            tx2gene <- tx2gene(
                organism,
                genomeBuild = genomeBuild,
                release = release
            )
        }
    }

    # Check for gene-to-symbol mappings
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = names(mcols(rowRanges)),
        severity = "warning"
    )

    # Generate rowData data.frame
    rowData <- as.data.frame(rowRanges)
    rownames(rowData) <- names(rowRanges)

    # Transcript to gene level counts (legacy) =================================
    # Now recommended to provide GTF in the bcbio run instead
    if (level == "transcript") {
        counts <- .transcriptToGeneLevelCounts(counts, tx2gene)
    }

    # Column data ==============================================================
    colData <- calculateMetrics(
        object = counts,
        rowData = rowData,
        prefilter = prefilter
    )
    # Bind the `nCount` column to the colData
    cbPass <- cbData[rownames(colData), "nCount", drop = FALSE]
    colData <- cbind(colData, cbPass)

    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the cells that passed prefiltering
        counts <- counts[, rownames(colData), drop = FALSE]
    }

    # Cell to sample mappings ==================================================
    cell2sample <- mapCellsToSamples(
        cells = rownames(colData),
        samples = rownames(sampleData)
    )

    # Metadata =================================================================
    metadata <- list(
        "version" = packageVersion,
        "pipeline" = pipeline,
        "uploadDir" = uploadDir,
        "sampleDirs" = sampleDirs,
        "sampleMetadataFile" = sampleMetadataFile,
        "sampleData" = sampleData,
        "interestingGroups" = interestingGroups,
        "cell2sample" = cell2sample,
        "organism" = organism,
        "genomeBuild" = genomeBuild,
        "ensemblRelease" = as.integer(ensemblRelease),
        "annotationHub" = as.list(ahMeta),
        "gffFile" = as.character(gffFile),
        "txdb" = txdb,
        "umiType" = umiType,
        "allSamples" = allSamples,
        "prefilter" = prefilter,
        # bcbio pipeline-specific ==============================================
        "projectDir" = projectDir,
        "template" = template,
        "runDate" = runDate,
        "yamlFile" = yamlFile,
        "yaml" = yaml,
        "tx2gene" = tx2gene,
        "dataVersions" = dataVersions,
        "programVersions" = programVersions,
        "bcbioLog" = bcbioLog,
        "bcbioCommandsLog" = bcbioCommandsLog,
        "cellularBarcodes" = cbList,
        "cellularBarcodeCutoff" = cellularBarcodeCutoff,
        "loadSingleCell" = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots) > 0L) {
        assert_are_disjoint_sets(names(metadata), names(dots))
        metadata <- c(metadata, dots)
    }

    # Return ===================================================================
    # Generate RangedSummarizedExperiment
    rse <- prepareSummarizedExperiment(
        assays = list(assay = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        isSpike = isSpike
    )
    # FIXME `as(rse, "SingleCellExperiment)` coercion doesn't return
    # objectVersion correctly in SingleCellExperiment v1.0.0
    sce <- SingleCellExperiment(
        assays = assays(rse),
        rowRanges = rowRanges(rse),
        colData = colData(rse),
        metadata = metadata(rse)
    )
    # Define spikeNames for spike-in sequences
    if (is.character(isSpike)) {
        for (i in seq_along(isSpike)) {
            isSpike(sce, isSpike[[i]]) <- isSpike[[i]]
        }
    }
    new("bcbioSingleCell", sce)
}
