#' @inherit bcbioSingleCell-class title description
#' @author Michael Steinbaugh
#' @note Updated 2019-08-27.
#' @export
#'
#' @inheritParams basejump::makeSingleCellExperiment
#' @inheritParams BiocParallel::bplapply
#' @inheritParams acidroxygen::params
#'
#' @section Remote data:
#'
#' When working in RStudio, we recommend connecting to the bcbio-nextgen run
#' directory as a remote connection over
#' [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @return `bcbioSingleCell`.
#'
#' @seealso
#' - `SingleCellExperiment::SingleCellExperiment()`.
#' - `.S4methods(class = "bcbioSingleCell")`.
#'
#' @examples
#' uploadDir <- system.file("extdata/indrops", package = "bcbioSingleCell")
#'
#' x <- bcbioSingleCell(uploadDir)
#' print(x)
#'
#' x <- bcbioSingleCell(
#'     uploadDir = uploadDir,
#'     sampleMetadataFile = file.path(uploadDir, "metadata.csv")
#' )
#' print(x)
bcbioSingleCell <- function(
    uploadDir,
    sampleMetadataFile = NULL,
    organism = NULL,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    gffFile = NULL,
    transgeneNames = NULL,
    spikeNames = NULL,
    interestingGroups = "sampleName",
    BPPARAM = BiocParallel::SerialParam(),  # nolint
    ...
) {
    sampleData <- NULL

    ## Legacy arguments --------------------------------------------------------
    ## nocov start
    dots <- list(...)
    call <- match.call()
    ## ensemblVersion
    if ("ensemblVersion" %in% names(call)) {
        stop("Use 'ensemblRelease' instead of 'ensemblVersion'.")
    }
    ## gtfFile
    if ("gtfFile" %in% names(call)) {
        stop("Use 'gffFile' instead of 'gtfFile'.")
    }
    ## annotable
    if ("annotable" %in% names(call)) {
        stop("Use 'gffFile' instead of 'annotable'.")
    }
    ## Error on unsupported arguments.
    assert(isSubset(
        x = setdiff(names(call), ""),
        y = names(formals())
    ))
    rm(call)
    ## nocov end

    ## Assert checks -----------------------------------------------------------
    assert(
        isADirectory(uploadDir),
        isString(sampleMetadataFile, nullOK = TRUE),
        isString(organism, nullOK = TRUE),
        isInt(ensemblRelease, nullOK = TRUE),
        isString(genomeBuild, nullOK = TRUE),
        isString(gffFile, nullOK = TRUE),
        isCharacter(transgeneNames, nullOK = TRUE),
        isCharacter(spikeNames, nullOK = TRUE),
        isCharacter(interestingGroups),
        identical(attr(class(BPPARAM), "package"), "BiocParallel")
    )
    if (isString(gffFile)) {
        isAFile(gffFile) || isAURL(gffFile)
    }

    ## Directory paths ---------------------------------------------------------
    uploadDir <- realpath(uploadDir)
    projectDir <- projectDir(uploadDir)
    sampleDirs <- sampleDirs(uploadDir)

    ## Sequencing lanes --------------------------------------------------------
    lanes <- detectLanes(sampleDirs)

    ## Project summary YAML ----------------------------------------------------
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- import(yamlFile)

    ## bcbio run information ---------------------------------------------------
    dataVersions <-
        importDataVersions(file.path(projectDir, "data_versions.csv"))
    assert(is(dataVersions, "DataFrame"))
    programVersions <-
        importProgramVersions(file.path(projectDir, "programs.txt"))
    assert(is(dataVersions, "DataFrame"))
    log <- import(file.path(projectDir, "bcbio-nextgen.log"))
    ## This step enables our minimal dataset inside the package to pass checks.
    tryCatch(
        expr = assert(isCharacter(log)),
        error = function(e) {
            message("'bcbio-nextgen.log' file is empty.")
        }
    )
    commandsLog <- import(file.path(projectDir, "bcbio-nextgen-commands.log"))
    ## This step enables our minimal dataset inside the package to pass checks.
    tryCatch(
        expr = assert(isCharacter(commandsLog)),
        error = function(e) {
            message("'bcbio-nextgen-commands.log' file is empty.")
        }
    )
    cutoff <- getBarcodeCutoffFromCommands(commandsLog)
    level <- getLevelFromCommands(commandsLog)
    umiType <- getUMITypeFromCommands(commandsLog)
    ## Check to see if we're dealing with a multiplexed platform.
    multiplexed <- any(vapply(
        X = c("dropseq", "indrop"),
        FUN = function(pattern) {
            grepl(pattern = pattern, x = umiType)
        },
        FUN.VALUE = logical(1L)
    ))

    ## User-defined sample metadata --------------------------------------------
    allSamples <- TRUE
    sampleData <- NULL
    if (isString(sampleMetadataFile)) {
        sampleData <- importSampleData(sampleMetadataFile, lanes = lanes)
        ## Error on incorrect reverse complement input.
        if ("sequence" %in% colnames(sampleData)) {
            sampleDirSequence <- str_extract(names(sampleDirs), "[ACGT]+$")
            if (identical(
                sort(sampleDirSequence),
                sort(as.character(sampleData[["sequence"]]))
            )) {
                stop(
                    "It appears that the reverse complement sequence of the ",
                    "i5 index barcodes were input into the sample metadata ",
                    "'sequence' column. bcbio outputs the revcomp into the ",
                    "sample directories, but the forward sequence should be ",
                    "used in the R package."
                )
            }
        }
        ## Allow sample selection by with this file.
        if (nrow(sampleData) < length(sampleDirs)) {
            sampleDirs <- sampleDirs[rownames(sampleData)]
            message(sprintf(
                fmt = "Loading a subset of samples: %s.",
                toString(basename(sampleDirs), width = 100L)
            ))
            allSamples <- FALSE
        }
    }

    ## Unfiltered cellular barcode distributions -------------------------------
    cbList <- .importReads(
        sampleDirs = sampleDirs,
        BPPARAM = BPPARAM
    )

    ## Assays ------------------------------------------------------------------
    ## Note that we're now allowing transcript-level counts.
    counts <- .importCounts(
        sampleDirs = sampleDirs,
        BPPARAM = BPPARAM
    )
    assert(hasValidDimnames(counts))

    ## Row data ----------------------------------------------------------------
    ## Annotation priority:
    ## 1. AnnotationHub.
    ##    - Requires `organism` to be declared.
    ##    - Ensure that Ensembl release and genome build match.
    ## 2. GTF/GFF file. Use the bcbio GTF if possible.
    ## 3. Fall back to slotting empty ranges. This is offered as support for
    ##    complex datasets (e.g. multiple organisms).
    if (isString(organism) && is.numeric(ensemblRelease)) {
        ## AnnotationHub (ensembldb).
        message("Using 'makeGRangesFromEnsembl()' for annotations.")
        rowRanges <- makeGRangesFromEnsembl(
            organism = organism,
            level = level,
            genomeBuild = genomeBuild,
            release = ensemblRelease
        )
    } else {
        ## GTF/GFF file.
        if (is.null(gffFile)) {
            ## Attempt to use bcbio GTF automatically.
            gffFile <- getGTFFileFromYAML(yaml)
        }
        if (!is.null(gffFile)) {
            message("Using 'makeGRangesFromGFF()' for annotations.")
            gffFile <- realpath(gffFile)
            rowRanges <- makeGRangesFromGFF(file = gffFile, level = level)
        } else {
            message("Slotting empty ranges into 'rowRanges()'.")
            rowRanges <- emptyRanges(rownames(counts))
        }
    }
    assert(is(rowRanges, "GRanges"))
    ## Attempt to get genome build and Ensembl release if not declared.
    ## Note that these will remain NULL when using GTF file (see above).
    if (is.null(genomeBuild)) {
        genomeBuild <- metadata(rowRanges)[["genomeBuild"]]
    }
    if (is.null(ensemblRelease)) {
        ensemblRelease <- metadata(rowRanges)[["ensemblRelease"]]
    }

    ## Column data -------------------------------------------------------------
    colData <- DataFrame(row.names = colnames(counts))
    ## Generate automatic sample metadata, if necessary.
    if (is.null(sampleData)) {
        if (isTRUE(multiplexed)) {
            ## Multiplexed samples without user-defined metadata.
            message(sprintf(
                fmt = paste0(
                    "'sampleMetadataFile' is recommended for ",
                    "multiplexed samples (e.g. %s)."
                ),
                umiType
            ))
            sampleData <- minimalSampleData(basename(sampleDirs))
        } else {
            sampleData <- getSampleDataFromYAML(yaml)
        }
    }
    assert(isSubset(rownames(sampleData), names(sampleDirs)))
    ## Join `sampleData` into cell-level `colData`.
    if (identical(nrow(sampleData), 1L)) {
        colData[["sampleID"]] <- as.factor(rownames(sampleData))
    } else {
        colData[["sampleID"]] <- mapCellsToSamples(
            cells = rownames(colData),
            samples = rownames(sampleData)
        )
    }
    sampleData[["sampleID"]] <- as.factor(rownames(sampleData))
    ## Need to ensure the `sampleID` factor levels match up, otherwise we'll get
    ## a warning during the `leftJjoin()` call below.
    assert(areSetEqual(
        x = levels(colData[["sampleID"]]),
        y = levels(sampleData[["sampleID"]])
    ))
    levels(sampleData[["sampleID"]]) <- levels(colData[["sampleID"]])
    colData <- leftJoin(colData, sampleData, by = "sampleID")
    assert(
        is(colData, "DataFrame"),
        hasRownames(colData)
    )

    ## Metadata ----------------------------------------------------------------
    runDate <- runDate(projectDir)
    interestingGroups <- camelCase(interestingGroups)
    assert(isSubset(interestingGroups, colnames(sampleData)))
    metadata <- list(
        allSamples = allSamples,
        bcbioCommandsLog = commandsLog,
        bcbioLog = log,
        call = standardizeCall(),
        cellularBarcodeCutoff = cutoff,
        cellularBarcodes = cbList,
        dataVersions = dataVersions,
        ensemblRelease = as.integer(ensemblRelease),
        genomeBuild = as.character(genomeBuild),
        gffFile = as.character(gffFile),
        interestingGroups = interestingGroups,
        lanes = lanes,
        level = level,
        organism = as.character(organism),
        pipeline = "bcbio",
        programVersions = programVersions,
        projectDir = projectDir,
        runDate = runDate,
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        umiType = umiType,
        uploadDir = uploadDir,
        version = .version,
        yaml = yaml
    )

    ## SingleCellExperiment ----------------------------------------------------
    object <- makeSingleCellExperiment(
        assays = SimpleList(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )

    ## Return ------------------------------------------------------------------
    ## Always prefilter, removing very low quality cells and/or genes.
    object <- calculateMetrics(object = object, prefilter = TRUE)
    ## Bind the `nRead` column into the cell metrics. These are the number of
    ## raw read counts prior to UMI disambiguation that bcbio uses for initial
    ## filtering (minimum_barcode_depth in YAML).
    colData <- colData(object)
    nRead <- .nRead(cbList)
    assert(
        is.integer(nRead),
        isSubset(rownames(colData), names(nRead)),
        areDisjointSets("nRead", colnames(colData))
    )
    ## Switched to "nRead" from "nCount" in v0.3.19.
    colData[["nRead"]] <- unname(nRead[rownames(colData)])
    colData <- colData[, sort(colnames(colData)), drop = FALSE]
    colData(object) <- colData
    new(Class = "bcbioSingleCell", object)
}
