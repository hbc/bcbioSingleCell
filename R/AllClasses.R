setOldClass(Classes = c("grouped_df", "tbl_df", "tibble"))



#' @rdname bcbioSingleCell
#' @aliases NULL
#' @exportClass bcbioSingleCell
#' @usage NULL
bcbioSingleCell <- setClass(
    Class = "bcbioSingleCell",
    contains = "SingleCellExperiment"
)



# Constructors =================================================================
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
#' @family S4 Class Definition
#' @docType class
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams bcbioBase::prepareSummarizedExperiment
#' @inheritParams general
#' @param uploadDir Path to final upload directory. This path is set when
#'   running "`bcbio_nextgen -w template`".
#' @param organism Organism name. Use the full latin name (e.g.
#'   "Homo sapiens"), since this will be input downstream to
#'   AnnotationHub and ensembldb, unless `gffFile` is set.
#' @param sampleMetadataFile Sample barcode metadata file. Optional for runs
#'   with demultiplixed index barcodes (e.g. SureCell), but otherwise required
#'   for runs with multipliexed FASTQs containing multiple index barcodes (e.g.
#'   inDrop).
#' @param ensemblRelease *Optional.* Ensembl release version. If `NULL`,
#'   defaults to current release, and does not typically need to be
#'   user-defined. Passed to AnnotationHub for `EnsDb` annotation matching,
#'   unless `gffFile` is set.
#' @param genomeBuild *Optional.* Ensembl genome build name (e.g. "GRCh38").
#'   This will be passed to AnnotationHub for `EnsDb` annotation matching,
#'   unless `gffFile` is set.
#' @param gffFile *Optional, not recommended.* By default, we recommend leaving
#'   this `NULL` for genomes that are supported on Ensembl. In this case, the
#'   row annotations ([rowRanges()]) will be obtained automatically from Ensembl
#'   by passing the `organism`, `genomeBuild`, and `ensemblRelease` arguments to
#'   AnnotationHub and ensembldb. For a genome that is not supported on Ensembl
#'   and/or AnnotationHub, a GFF/GTF (General Feature Format) file is required.
#'   Generally, we recommend using a GTF (GFFv2) file here over a GFF3 file if
#'   possible, although all GFF formats are supported. The function will
#'   internally generate a `TxDb` containing transcript-to-gene mappings and
#'   construct a `GRanges` object containing the genomic ranges ([rowRanges()]).
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
    organism,
    sampleMetadataFile,
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

    # Legacy arguments =========================================================
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
    # organism missing
    if (!"organism" %in% names(call)) {
        stop("`organism` is now required")
    }
    # gtfFile
    if ("gtfFile" %in% names(call)) {
        warning("Use `gffFile` instead of `gtfFile`")
        gffFile <- call[["gtfFile"]]
        dots[["gtfFile"]] <- NULL
    }
    dots <- Filter(Negate(is.null), dots)

    # Assert checks ============================================================
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    if (missing(sampleMetadataFile)) {
        sampleMetadataFile <- NULL
    }
    assertIsAStringOrNULL(sampleMetadataFile)
    assert_is_character(interestingGroups)
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsAStringOrNULL(genomeBuild)
    assertIsCharacterOrNULL(transgeneNames)
    assertIsCharacterOrNULL(spikeNames)
    assertIsAStringOrNULL(gffFile)
    if (is_a_string(gffFile)) {
        assert_all_are_existing_files(gffFile)
    }

    # Directory paths ==========================================================
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    projectDir <- dir(
        uploadDir,
        pattern = bcbioBase::projectDirPattern,
        full.names = FALSE,
        recursive = FALSE
    )
    assert_is_a_string(projectDir)
    message(projectDir)
    match <- str_match(projectDir, bcbioBase::projectDirPattern)
    runDate <- as.Date(match[[2L]])
    template <- match[[3L]]
    projectDir <- file.path(uploadDir, projectDir)
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sequencing lanes =========================================================
    if (any(grepl(bcbioBase::lanePattern, sampleDirs))) {
        lanes <- str_match(names(sampleDirs), bcbioBase::lanePattern) %>%
            .[, 2L] %>%
            unique() %>%
            length()
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"
        ))
    } else {
        lanes <- 1L
    }

    # Project summary YAML =====================================================
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- readYAML(yamlFile)

    # bcbio run information ====================================================
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

    # Cellular barcode cutoff ==================================================
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

    # Detect gene or transcript-level output ===================================
    genemapPattern <- "--genemap (.+)-tx2gene.tsv"
    if (any(grepl(genemapPattern, bcbioCommandsLog))) {
        level <- "genes"
    } else {
        level <- "transcripts"
    }

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
    message(paste("UMI type:", umiType))

    # Sample metadata ==========================================================
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

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    assertFormalInterestingGroups(sampleData, interestingGroups)

    # Subset sample directories by metadata ====================================
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

    # Assays ===================================================================
    message(paste("Reading counts as", level))
    countsList <- .sparseCountsList(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        umiType = umiType
    )
    # Ensure samples with empty matrices (`NULL`) are filtered
    countsList <- Filter(Negate(is.null), countsList)
    counts <- do.call(cbind, countsList)

    # Require transcript to gene conversion (legacy) ===========================
    if (level == "transcripts") {
        message("Converting transcripts to genes")

        if (!is_a_string(gffFile)) {
            stop("GFF is required to convert transcripts to genes")
        }

        tx2gene <- makeTx2geneFromGFF(gffFile)

        # Add spike-ins to tx2gene, if necessary
        if (is.character(isSpike)) {
            assert_are_disjoint_sets(rownames(tx2gene), isSpike)
            spike <- data.frame(
                "txID" = isSpike,
                "geneID" = isSpike,
                row.names = isSpike,
                stringsAsFactors = FALSE
            )
            tx2gene <- rbind(spike, tx2gene)
        }

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

    # Unfiltered cellular barcode distributions ================================
    cbList <- .cellularBarcodesList(sampleDirs)
    cbData <- .bindCellularBarcodes(cbList)

    # Row data =================================================================
    rowRangesMetadata <- NULL
    if (is_a_string(gffFile)) {
        rowRanges <- makeGRangesFromGFF(gffFile, format = "genes")
    } else if (is_a_string(organism)) {
        # ah: AnnotationHub
        ah <- makeGRangesFromEnsembl(
            organism = organism,
            format = "genes",
            genomeBuild = genomeBuild,
            release = ensemblRelease,
            metadata = TRUE
        )
        assert_is_list(ah)
        assert_are_identical(names(ah), c("data", "metadata"))
        rowRanges <- ah[["data"]]
        assert_is_all_of(rowRanges, "GRanges")
        rowRangesMetadata <- ah[["metadata"]]
        assert_is_data.frame(rowRangesMetadata)
    } else {
        rowRanges <- emptyRanges(rownames(counts))
    }

    rowData <- as.data.frame(rowRanges)
    rownames(rowData) <- names(rowRanges)

    # Column data ==============================================================
    # Always prefilter, removing very low quality cells with no UMIs or genes
    metrics <- metrics(counts, rowData = rowData, prefilter = TRUE)
    colData <- as(metrics, "DataFrame")

    # Cell to sample mappings
    cell2sample <- mapCellsToSamples(
        cells = rownames(colData),
        samples = rownames(sampleData)
    )

    colData[["cellID"]] <- rownames(colData)
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

    # Subset the counts to match the prefiltered metrics
    counts <- counts[, rownames(colData), drop = FALSE]

    # Metadata =================================================================
    metadata <- list(
        "version" = packageVersion,
        "pipeline" = pipeline,
        "level" = level,
        "uploadDir" = uploadDir,
        "sampleDirs" = sampleDirs,
        "sampleMetadataFile" = as.character(sampleMetadataFile),
        "interestingGroups" = interestingGroups,
        "organism" = as.character(organism),
        "genomeBuild" = as.character(genomeBuild),
        "ensemblRelease" = as.integer(ensemblRelease),
        "rowRangesMetadata" = rowRangesMetadata,
        "sampleData" = sampleData,
        "cell2sample" = as.factor(cell2sample),
        "umiType" = umiType,
        "allSamples" = allSamples,
        # bcbio pipeline-specific ----------------------------------------------
        "projectDir" = projectDir,
        "template" = template,
        "runDate" = runDate,
        "yaml" = yaml,
        "gffFile" = as.character(gffFile),
        "tx2gene" = tx2gene,
        "dataVersions" = dataVersions,
        "programVersions" = programVersions,
        "bcbioLog" = bcbioLog,
        "bcbioCommandsLog" = bcbioCommandsLog,
        "cellularBarcodes" = cbList,
        "cellularBarcodeCutoff" = cellularBarcodeCutoff,
        "call" = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots)) {
        assert_are_disjoint_sets(names(metadata), names(dots))
        metadata <- c(metadata, dots)
    }

    # Return ===================================================================
    .new.bcbioSingleCell(
        assays = list("counts" = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}



# Validity Checks ==============================================================
setValidity(
    "bcbioSingleCell",
    function(object) {
        stopifnot(metadata(object)[["version"]] >= 0.1)
        stopifnot(!.hasSlot(object, "bcbio"))

        # Assays ===============================================================
        assert_are_identical("counts", names(assays(object)))

        # Row data =============================================================
        assert_is_all_of(rowRanges(object), "GRanges")
        assert_is_all_of(rowData(object), "DataFrame")

        # Column data ==========================================================
        # TODO Inform the user about `colData` structure change made in v0.1.7
        # assert_is_subset(
        #     x = colnames(metadata(object)[["sampleData"]]),
        #     y = colnames(colData(object))
        # )

        # Metadata =============================================================
        metadata <- metadata(object)

        # Optional metadata:
        # - filterCells
        # - filterGenes
        # - filterParams
        # - filterSummary
        # - lanes: integer
        # - rowRangesMetadata: tbl_df
        # - tx2gene: data.frame
        #
        # bcbio-specific:
        # - bcbioCommandsLog: character
        # - bcbioLog: character
        # - dataVersions: tbl_df
        # - gffFile: character
        # - programVersions: tbl_df
        # - projectDir: character
        # - runDate: Date
        # - template: character
        # - yaml: list

        # Class checks
        requiredMetadata <- list(
            "allSamples" = "logical",
            "cell2sample" = "factor",
            "date" = "Date",
            "devtoolsSessionInfo" = "session_info",
            "ensemblRelease" = "integer",
            "genomeBuild" = "character",
            "interestingGroups" = "character",
            "level" = "character",
            "organism" = "character",
            "pipeline" = "character",
            "sampleData" = "data.frame",
            "sampleDirs" = "character",
            "sampleMetadataFile" = "character",
            "umiType" = "character",
            "uploadDir" = "character",
            "utilsSessionInfo" = "sessionInfo",
            "version" = "package_version",
            "wd" = "character"
        )
        classChecks <- invisible(mapply(
            name <- names(requiredMetadata),
            expected <- requiredMetadata,
            MoreArgs = list(metadata = metadata),
            FUN = function(name, expected, metadata) {
                actual <- class(metadata[[name]])
                if (!length(intersect(expected, actual))) {
                    FALSE
                } else {
                    TRUE
                }
            },
            SIMPLIFY = TRUE,
            USE.NAMES = TRUE
        ))
        if (!all(classChecks)) {
            print(classChecks)
            stop(paste(
                "Metadata class checks failed.",
                bcbioBase::updateMessage,
                sep = "\n"
            ))
        }

        # level
        assert_is_subset(
            x = metadata[["level"]],
            y = c("genes", "transcripts")
        )

        TRUE
    }
)
