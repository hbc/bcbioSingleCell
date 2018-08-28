#' `bcbioSingleCell` Object and Constructor
#'
#' `bcbioSingleCell` is an S4 class that extends `SingleCellExperiment`, and is
#' designed to store a bcbio single-cell RNA-seq analysis. This class contains
#' read counts saved as a sparse matrix (`dgCMatrix`), sample metadata, and cell
#' quality control metrics.
#'
#' @section Remote Data:
#'
#' When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @note `bcbioSingleCell` extended `SummarizedExperiment` prior to v0.1.0,
#'   where we migrated to `SingleCellExperiment`.
#'
#' @rdname bcbioSingleCell
#' @aliases bcbioSingleCell-class
#' @family S4 Object
#' @docType class
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams basejump::makeSummarizedExperiment
#' @inheritParams general
#' @param uploadDir `string`. Path to final upload directory. This path is set
#'   when running "`bcbio_nextgen -w template`".
#' @param organism `string` or `NULL`. Organism name. Use the full Latin name
#'   (e.g. "Homo sapiens"), since this will be input downstream to AnnotationHub
#'   and ensembldb, unless `gffFile` is set. If left `NULL` (*not recommended*),
#'   the function call will skip loading gene-level annotations into
#'   [rowRanges()]. This can be useful for poorly annotation genomes or
#'   experiments involving multiple genomes.
#' @param sampleMetadataFile `string` or `NULL`. Sample barcode metadata file.
#'   Optional for runs with demultiplixed index barcodes (e.g. SureCell) and can
#'   be left `NULL`, but otherwise required for runs with multipliexed FASTQs
#'   containing multiple index barcodes (e.g. inDrop).
#' @param ensemblRelease `scalar integer` or `NULL`. Ensembl release version. If
#'   `NULL`, defaults to current release, and does not typically need to be
#'   user-defined. Passed to AnnotationHub for `EnsDb` annotation matching,
#'   unless `gffFile` is set.
#' @param genomeBuild `string` or `NULL`. Ensembl genome build name (e.g.
#'   "GRCh38"). This will be passed to AnnotationHub for `EnsDb` annotation
#'   matching, unless `gffFile` is set.
#' @param gffFile `string` or `NULL`. *Advanced use only; not recommended.* By
#'   default, we recommend leaving this `NULL` for genomes that are supported on
#'   Ensembl. In this case, the row annotations ([rowRanges()]) will be obtained
#'   automatically from Ensembl by passing the `organism`, `genomeBuild`, and
#'   `ensemblRelease` arguments to AnnotationHub and ensembldb. For a genome
#'   that is not supported on Ensembl and/or AnnotationHub, a GFF/GTF (General
#'   Feature Format) file is required. Generally, we recommend using a GTF
#'   (GFFv2) file here over a GFF3 file if possible, although all GFF formats
#'   are supported. The function will internally generate transcript-to-gene
#'   mappings and construct a `GRanges` object containing the genomic ranges
#'   ([rowRanges()]).
#' @param ... Additional arguments, to be stashed in the [metadata()] slot.
#'
#' @return `bcbioSingleCell`.
#' @export
#'
#' @seealso
#' - [SingleCellExperiment::SingleCellExperiment()].
#' - `.S4methods(class = "bcbioSingleCell")`.
#'
#' @examples
#' uploadDir <- system.file("extdata/indrops", package = "bcbioSingleCell")
#' x <- bcbioSingleCell(
#'     uploadDir = uploadDir,
#'     organism = "Homo sapiens",
#'     sampleMetadataFile = file.path(uploadDir, "metadata.csv"),
#'     ensemblRelease = 87L
#' )
#' show(x)
#' validObject(x)
bcbioSingleCell <- function(
    uploadDir,
    organism = NULL,
    sampleMetadataFile = NULL,
    interestingGroups = "sampleName",
    ensemblRelease = NULL,
    genomeBuild = NULL,
    transgeneNames = NULL,
    spikeNames = NULL,
    gffFile = NULL,
    ...
) {
    dots <- list(...)
    pipeline <- "bcbio"

    # Legacy arguments ---------------------------------------------------------
    call <- match.call(expand.dots = TRUE)
    # annotable
    if ("annotable" %in% names(call)) {
        stop("Use `gffFile` instead of `annotable`")
    }
    # ensemblVersion
    if ("ensemblVersion" %in% names(call)) {
        warning("Use `ensemblRelease` instead of `ensemblVersion`")
        ensemblRelease <- call[["ensemblVersion"]]
        dots[["ensemblVersion"]] <- NULL
    }
    # gtfFile
    if ("gtfFile" %in% names(call)) {
        warning("Use `gffFile` instead of `gtfFile`")
        gffFile <- call[["gtfFile"]]
        dots[["gtfFile"]] <- NULL
    }
    # organism
    if (!"organism" %in% names(call)) {
        message("`organism` is now recommended, to acquire gene annotations")
    }
    dots <- Filter(Negate(is.null), dots)

    # Assert checks ------------------------------------------------------------
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    assertIsAStringOrNULL(sampleMetadataFile)
    assert_is_character(interestingGroups)
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsAStringOrNULL(genomeBuild)
    assert_is_any_of(transgeneNames, c("character", "NULL"))
    assert_is_any_of(spikeNames, c("character", "NULL"))
    assertIsAStringOrNULL(gffFile)
    if (is_a_string(gffFile)) {
        assert_all_are_existing_files(gffFile)
    }

    # Directory paths ----------------------------------------------------------
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    projectDir <- projectDir(uploadDir)
    sampleDirs <- sampleDirs(uploadDir)

    # Run date and template name -----------------------------------------------
    # Get run date and template name from project directory.
    # This information will be stashed in `metadata()`.
    match <- str_match(
        string = basename(projectDir),
        pattern = projectDirPattern
    )
    runDate <- as.Date(match[[2L]])
    template <- match[[3L]]
    rm(match)

    # Sequencing lanes ---------------------------------------------------------
    if (any(grepl(lanePattern, sampleDirs))) {
        lanes <- str_match(names(sampleDirs), lanePattern) %>%
            .[, 2L] %>%
            unique() %>%
            length()
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"
        ))
    } else {
        lanes <- 1L
    }

    # Project summary YAML -----------------------------------------------------
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- readYAML(yamlFile)

    # bcbio run information ----------------------------------------------------
    dataVersions <- readDataVersions(
        file = file.path(projectDir, "data_versions.csv")
    )
    assert_is_tbl_df(dataVersions)

    programVersions <- readProgramVersions(
        file = file.path(projectDir, "programs.txt")
    )
    assert_is_tbl_df(programVersions)

    bcbioLog <- readLog(
        file = file.path(projectDir, "bcbio-nextgen.log")
    )
    assert_is_character(bcbioLog)

    bcbioCommandsLog <- readLog(
        file = file.path(projectDir, "bcbio-nextgen-commands.log")
    )
    assert_is_character(bcbioCommandsLog)

    # Cellular barcode cutoff --------------------------------------------------
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

    message(paste(
        cellularBarcodeCutoff,
        "reads per cellular barcode cutoff detected"
    ))

    # Detect gene or transcript-level output -----------------------------------
    genemapPattern <- "--genemap (.+)-tx2gene.tsv"
    if (any(grepl(genemapPattern, bcbioCommandsLog))) {
        level <- "genes"
    } else {
        level <- "transcripts"
    }

    # Molecular barcode (UMI) type ---------------------------------------------
    umiPattern <- "fastqtransform.*/(.*)\\.json"
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
    message(paste("UMI type:", umiType))

    # Sample metadata ----------------------------------------------------------
    # External file required for inDrop
    if (grepl("indrop", umiType) && is.null(sampleMetadataFile)) {
        stop(paste(
            "inDrop samples require `sampleMetadataFile`",
            "containing the index barcode sequences"
        ))
    }

    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleData(sampleMetadataFile)
    } else {
        sampleData <- readYAMLSampleData(yamlFile)
    }

    # Check for incorrect reverse complement input
    if ("sequence" %in% colnames(sampleData)) {
        sampleDirSequence <- str_extract(names(sampleDirs), "[ACGT]+$")
        if (identical(
            sort(sampleDirSequence),
            sort(as.character(sampleData[["sequence"]]))
        )) {
            stop(paste(
                "It appears that the reverse complement sequence of the",
                "i5 index barcodes were input into the sample metadata",
                "`sequence` column. bcbio outputs the revcomp into the",
                "sample directories, but the forward sequence should be",
                "used in the R package."
            ))
        }
    }

    assert_is_subset(rownames(sampleData), names(sampleDirs))
    sampleData <- sanitizeSampleData(sampleData)

    # Interesting groups -------------------------------------------------------
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    assertFormalInterestingGroups(sampleData, interestingGroups)

    # Subset sample directories by metadata ------------------------------------
    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (nrow(sampleData) < length(sampleDirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs[rownames(sampleData)]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Assays -------------------------------------------------------------------
    message(paste("Reading counts as", level))
    counts <- .readCounts(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        format = "mtx"
    )

    # Require transcript to gene conversion (legacy) ---------------------------
    if (level == "transcripts") {
        message("Converting transcripts to genes")

        # GFF file is required
        if (!is_a_string(gffFile)) {
            stop("GFF is required to convert transcripts to genes")
        }

        # Note that this data.frame won't contain transcript versions
        tx2gene <- makeTx2geneFromGFF(gffFile)

        # Add spike-ins to tx2gene, if necessary
        if (is.character(isSpike)) {
            assert_are_disjoint_sets(rownames(tx2gene), isSpike)
            spike <- data.frame(
                transcriptID = isSpike,
                geneID = isSpike,
                row.names = isSpike,
                stringsAsFactors = FALSE
            )
            tx2gene <- rbind(spike, tx2gene)
        }

        # Ensure Ensembl transcript versions are removed
        counts <- stripTranscriptVersions(counts)

        # Resize the tx2gene to match the matrix rownames
        assert_is_subset(rownames(counts), rownames(tx2gene))
        tx2gene <- tx2gene[rownames(counts), , drop = FALSE]

        # Now we're ready to assign `geneID` and aggregate
        rownames(counts) <- tx2gene[["geneID"]]
        counts <- aggregate.Matrix(
            x = counts,
            groupings = rownames(counts),
            fun = "sum"
        )
        level <- "genes"
    } else {
        tx2gene <- NULL
    }

    # Unfiltered cellular barcode distributions --------------------------------
    cbList <- .cellularBarcodesList(sampleDirs)
    cbData <- .bindCellularBarcodes(cbList)

    # Row data -----------------------------------------------------------------
    rowRangesMetadata <- NULL
    if (is_a_string(gffFile)) {
        rowRanges <- makeGRangesFromGFF(gffFile, format = "genes")
    } else if (is_a_string(organism)) {
        # ah: AnnotationHub
        ah <- makeGRangesFromEnsembl(
            organism = organism,
            format = "genes",
            build = genomeBuild,
            release = ensemblRelease,
            metadata = TRUE
        )
        assert_is_list(ah)
        assert_are_identical(names(ah), c("data", "metadata"))
        rowRanges <- ah[["data"]]
        assert_is_all_of(rowRanges, "GRanges")
        rowRangesMetadata <- ah[["metadata"]]
        assert_is_data.frame(rowRangesMetadata)
        genomeBuild <- rowRangesMetadata %>%
            filter(!!sym("name") == "genome_build") %>%
            pull("value")
        assert_is_a_string(genomeBuild)
    } else {
        rowRanges <- emptyRanges(rownames(counts))
    }
    assert_is_subset(rownames(counts), names(rowRanges))

    # Column data --------------------------------------------------------------
    # Always prefilter, removing very low quality cells with no UMIs or genes
    metrics <- .metrics(
        object = counts,
        rowRanges = rowRanges,
        prefilter = TRUE
    )

    # Subset the counts to match the prefiltered metrics
    counts <- counts[, rownames(metrics), drop = FALSE]

    colData <- as(metrics, "DataFrame")
    colData[["cellID"]] <- rownames(colData)
    cell2sample <- .mapCellsToSamples(
        cells = rownames(colData),
        samples = rownames(sampleData)
    )
    colData[["sampleID"]] <- cell2sample
    sampleData[["sampleID"]] <- rownames(sampleData)
    colData <- merge(
        x = colData,
        y = sampleData,
        by = "sampleID",
        all.x = TRUE
    )
    rownames(colData) <- colData[["cellID"]]
    colData[["cellID"]] <- NULL
    sampleData[["sampleID"]] <- NULL

    # Bind the `nCount` column into the colData
    nCount <- cbData[rownames(colData), "nCount", drop = FALSE]
    colData <- cbind(nCount, colData)

    # Metadata -----------------------------------------------------------------
    metadata <- list(
        version = packageVersion,
        pipeline = pipeline,
        level = level,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        interestingGroups = interestingGroups,
        organism = as.character(organism),
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        rowRangesMetadata = rowRangesMetadata,
        sampleData = as(sampleData, "DataFrame"),
        cell2sample = as.factor(cell2sample),
        umiType = umiType,
        allSamples = allSamples,
        # bcbio pipeline-specific ----------------------------------------------
        projectDir = projectDir,
        template = template,
        runDate = runDate,
        yaml = yaml,
        gffFile = as.character(gffFile),
        tx2gene = tx2gene,
        dataVersions = dataVersions,
        programVersions = programVersions,
        bcbioLog = bcbioLog,
        bcbioCommandsLog = bcbioCommandsLog,
        cellularBarcodes = cbList,
        cellularBarcodeCutoff = cellularBarcodeCutoff,
        call = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots)) {
        assert_are_disjoint_sets(names(metadata), names(dots))
        metadata <- c(metadata, dots)
    }

    # Return -------------------------------------------------------------------
    .new.bcbioSingleCell(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}
