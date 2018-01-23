#' Load bcbio Single-Cell RNA-Seq Data
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom bcbioBase annotable camel detectOrganism gene2symbolFromGTF
#'   prepareSummarizedExperiment readDataVersions readGTF readLogFile
#'   readProgramVersions readSampleMetadataFile readYAML sampleYAMLMetadata
#'   tx2geneFromGTF
#' @importFrom Matrix cBind
#' @importFrom pbapply pblapply
#' @importFrom rlang is_string
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
#' @param organism *Optional*. Organism name. Use the full latin name (e.g.
#'   "Homo sapiens"), since this will be input downstream to
#'   AnnotationHub/ensembldb. If set, this genome must be supported on Ensembl.
#'   Normally this can be left `NULL`, and the function will attempt to detect
#'   the organism automatically using [detectOrganism()].
#' @param ensemblVersion *Optional*. Ensembl release version. If `NULL`,
#'   defaults to current release, and does not typically need to be
#'   user-defined. This parameter can be useful for matching Ensembl annotations
#'   against an outdated bcbio annotation build.
#' @param genomeBuild *Optional*. Genome build. Normally this can be left `NULL`
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
#' # Minimal working example dataset
#' extdataDir <- system.file("extdata", package = "bcbioSingleCell")
#' uploadDir <- file.path(extdataDir, "harvard_indrop_v3")
#' sampleMetadataFile <- file.path(extdataDir, "harvard_indrop_v3.xlsx")
#' bcb <- loadSingleCell(
#'     uploadDir = uploadDir,
#'     sampleMetadataFile = sampleMetadataFile)
#'
#' # Mus musculus
#' # Run with Ensembl 88 transcriptome FASTA and GTF files
#' \dontrun{
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
    sampleMetadataFile = NULL,
    interestingGroups = "sampleName",
    gtfFile = NULL,
    annotable = TRUE,
    organism = NULL,
    ensemblVersion = NULL,
    genomeBuild = NULL,
    prefilter = TRUE,
    ...) {
    pipeline <- "bcbio"

    # Parameter integrity checks ===============================================
    # uploadDir
    if (!is_string(uploadDir)) {
        abort("`uploadDir` must be string")
    }
    # sampleMetadataFile
    if (!any(
        is_string(sampleMetadataFile),
        is.null(sampleMetadataFile)
    )) {
        abort("`sampleMetadataFile` must be string or NULL")
    }
    # interestingGroups
    if (!is.character(interestingGroups)) {
        abort("`interestingGroups` must be character")
    }
    # gtfFile
    if (!any(
        is_string(gtfFile),
        is.null(gtfFile)
    )) {
        abort("`gtfFile` must be string or NULL")
    }
    # annotable
    if (!any(
        is.logical(annotable),
        is.data.frame(annotable),
        is.null(annotable)
    )) {
        abort("`annotable` must be logical, data.frame, or NULL")
    }
    # organism
    if (!any(
        is_string(organism),
        is.null(organism)
    )) {
        abort("`organism` must be string or NULL")
    }
    # ensemblVersion
    if (!any(
        is.null(ensemblVersion),
        is.numeric(ensemblVersion) && length(ensemblVersion) == 1
    )) {
        abort("`ensemblVersion` must be single numeric or NULL")
    }
    # genomeBuild
    if (!any(
        is_string(genomeBuild),
        is.null(genomeBuild)
    )) {
        abort("`genomeBuild` must be string or NULL")
    }
    # prefilter
    if (!is.logical(prefilter)) {
        abort("`prefilter` must be logical")
    }

    # Directory paths ==========================================================
    # Check connection to final upload directory
    if (!dir.exists(uploadDir)) {
        abort("Final upload directory does not exist")
    }
    uploadDir <- normalizePath(uploadDir)
    projectDir <- dir(
        uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE)
    if (length(projectDir) != 1) {
        abort("Failed to detect project directory")
    }
    inform(projectDir)
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
        inform(paste(
            lanes, "sequencing lane detected", "(technical replicates)"))
    } else {
        lanes <- 1
    }

    # Project summary YAML =====================================================
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- readYAML(yamlFile)

    # Log files ================================================================
    inform("Reading log files")
    bcbioLog <- readLogFile(
        file.path(projectDir, "bcbio-nextgen.log"))
    bcbioCommandsLog <- readLogFile(
        file.path(projectDir, "bcbio-nextgen-commands.log"))

    # Cellular barcode cutoff
    cellularBarcodeCutoffPattern <- "--cb_cutoff (\\d+)"
    if (length(bcbioCommandsLog)) {
        match <- str_match(
            string = bcbioCommandsLog,
            pattern = cellularBarcodeCutoffPattern)
        cellularBarcodeCutoff <- match %>%
            .[, 2L] %>%
            na.omit() %>%
            unique() %>%
            as.numeric()
    } else {
        cellularBarcodeCutoff <- NULL
    }
    if (!is.null(cellularBarcodeCutoff)) {
        inform(paste(
            cellularBarcodeCutoff,
            "reads per cellular barcode cutoff detected"
        ))
    }

    # Detect MatrixMarket output at transcript or gene level
    if (is.character(bcbioCommandsLog)) {
        # This grep pattern may not be strict enough against the file path
        genemapPattern <- "--genemap (.+)-tx2gene.tsv"
        if (any(grepl(x = bcbioCommandsLog, pattern = genemapPattern))) {
            countsLevel <- "gene"
        } else {
            countsLevel <- "transcript"
        }
    } else {
        # The pipeline now defaults to gene level. If there's no log file,
        # let's assume the counts are genes.
        countsLevel <- "gene"
    }

    # Data versions and programs ===============================================
    inform("Reading data and program versions")
    dataVersions <- readDataVersions(
        file.path(projectDir, "data_versions.csv"))
    programs <- readProgramVersions(
        file.path(projectDir, "programs.txt"))

    # Detect genome build
    if (!is_string(genomeBuild)) {
        if (is.data.frame(dataVersions)) {
            genomeBuild <- dataVersions %>%
                .[.[["resource"]] == "transcripts", "genome", drop = TRUE]
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
                warn("Genome detection from bcbio commands failed")
                genomeBuild <- NULL
            }
        }
    }

    # Detect organism
    if (!is_string(organism)) {
        if (!is_string(genomeBuild)) {
            abort("Organism detection by genome build failed")
        }
        organism <- detectOrganism(genomeBuild)
    }
    if (!is_string(organism)) {
        abort("Invalid organism")
    }
    inform(paste(
        paste("Organism:", organism),
        paste("Genome build:", genomeBuild),
        sep = "\n"
    ))

    # Molecular barcode (UMI) type =============================================
    if (is.character(bcbioCommandsLog)) {
        umiPattern <- "/umis/([a-z0-9\\-]+)\\.json"
        if (any(grepl(x = bcbioCommandsLog, pattern = umiPattern))) {
            umiType <- str_match(bcbioCommandsLog, umiPattern) %>%
                .[, 2] %>%
                na.omit() %>%
                unique() %>%
                gsub(x = .,
                     pattern = "-transform",
                     replacement = "")
            inform(paste("UMI type:", umiType))
        } else {
            warn(paste(
                "Failed to detect UMI type from commands log JSON file grep"
            ))
            umiType <- NULL
        }
    } else {
        umiType <- NULL
    }
    # Assume samples are inDrop platform, by default
    if (is.null(umiType)) {
        warn(paste(
            "Assuming unknown UMI is 'harvard-indrop-v3' (default)"
        ))
        umiType <- "harvard-indrop-v3"
    }

    # Sample metadata ==========================================================
    if (is_string(sampleMetadataFile)) {
        if (!file.exists(sampleMetadataFile)) {
            abort("'sampleMetadataFile missing")
        }
        sampleMetadataFile <- normalizePath(sampleMetadataFile)
        sampleMetadata <- readSampleMetadataFile(sampleMetadataFile)
    } else {
        if (grepl(x = umiType, pattern = "indrop")) {
            # Enforce `sampleMetadataFile` for multiplexed data containing
            # index barcodes (e.g. inDrop)
            abort("`sampleMetadataFile` is required for inDrop samples")
        }
        sampleMetadata <- sampleYAMLMetadata(yaml)
    }
    # Check that `sampleID` matches `sampleDirs`
    if (!all(rownames(sampleMetadata) %in% names(sampleDirs))) {
        # Check for flipped index sequence. Sample metadata should use the
        # forward sequence (`sequence`), whereas bcbio names the sample
        # directories with the reverse complement (`revcomp`).
        if ("sequence" %in% colnames(sampleMetadata)) {
            sampleDirSequence <- str_extract(
                string = names(sampleDirs),
                pattern = "[ACGT]+$")
            if (identical(
                sort(sampleDirSequence),
                sort(as.character(sampleMetadata[["sequence"]]))
            )) {
                warn(paste(
                    "It appears that the reverse complement sequence of the",
                    "i5 index barcode(s) was input into the sample metadata",
                    "'sequence' column. bcbio outputs the revcomp into the",
                    "sample directories, but the forward sequence should be",
                    "used in the R package."
                ))
            }
        }
        abort("Sample directory names don't match the sample metadata file")
    }

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    # Check to ensure interesting groups are defined
    if (!all(interestingGroups %in% colnames(sampleMetadata))) {
        abort("Interesting groups missing in sample metadata")
    }

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
        annotable <- annotable(annotable, uniqueSymbol = FALSE)
    } else {
        warn("Loading run without gene annotable")
        annotable <- NULL
    }

    # GTF annotations
    if (is_string(gtfFile)) {
        gtfFile <- normalizePath(gtfFile)
        gtf <- readGTF(gtfFile)
    } else {
        gtf <- NULL
    }

    # Transcript-to-gene mappings
    if (countsLevel == "transcript") {
        if (!is.data.frame(gtf)) {
            abort(paste(
                "GTF required to convert transcript-level counts to gene-level"
            ))
        }
        tx2gene <- tx2geneFromGTF(gtf)
    } else {
        # Not applicable to gene-level bcbio output
        tx2gene <- NA
    }

    # Gene-to-symbol mappings
    if (is_string(gtfFile)) {
        gene2symbol <- gene2symbolFromGTF(gtf)
    } else if (is.data.frame(annotable)) {
        gene2symbol <- annotable[, c("ensgene", "symbol")]
    } else {
        warn("Loading run without gene-to-symbol mappings")
        gene2symbol <- NULL
    }

    # Cellular barcode distributions ===========================================
    cbList <- .cellularBarcodesList(sampleDirs)
    cbData <- .bindCellularBarcodes(cbList)

    # Counts ===================================================================
    inform(paste("Reading counts at", countsLevel, "level"))
    # Migrate this to `mapply()` method in future update
    sparseList <- pblapply(seq_along(sampleDirs), function(a) {
        .readSparseCounts(
            sampleDirs[a],
            pipeline = pipeline,
            umiType = umiType)
    })
    names(sparseList) <- names(sampleDirs)
    # Combine the individual per-sample transcript-level sparse matrices into a
    # single sparse matrix
    counts <- do.call(Matrix::cBind, sparseList)
    # Convert counts from transcript-level to gene-level, if necessary
    if (countsLevel == "transcript") {
        counts <- .transcriptToGeneLevelCounts(counts, tx2gene)
    }

    # Metrics ==================================================================
    metrics <- calculateMetrics(
        counts,
        annotable = annotable,
        prefilter = prefilter)
    # Bind the `nCount` column to the metrics
    cbPass <- cbData[rownames(metrics), "nCount", drop = FALSE]
    metrics <- cbind(metrics, cbPass)

    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the cells that passed prefiltering
        counts <- counts[, rownames(metrics)]
    }

    # Cell to sample mappings ==================================================
    cell2sample <- cell2sample(
        rownames(metrics),
        samples = rownames(sampleMetadata)
    )

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
    if (length(dots) > 0) {
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
    new("bcbioSingleCell",
        se,
        bcbio = as(bcbio, "SimpleList"))
}
