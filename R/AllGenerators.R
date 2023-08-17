#' @inherit bcbioSingleCell-class title description
#' @author Michael Steinbaugh
#' @note Updated 2023-08-17.
#' @export
#'
#' @inheritParams AcidSingleCell::makeSingleCellExperiment
#' @inheritParams AcidRoxygen::params
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
bcbioSingleCell <-
    function(uploadDir,
             sampleMetadataFile = NULL,
             organism = NULL,
             ensemblRelease = NULL,
             genomeBuild = NULL,
             gffFile = NULL,
             transgeneNames = NULL,
             interestingGroups = "sampleName") {
        assert(
            isADirectory(uploadDir),
            isString(sampleMetadataFile, nullOK = TRUE),
            isString(organism, nullOK = TRUE),
            isInt(ensemblRelease, nullOK = TRUE),
            isString(genomeBuild, nullOK = TRUE),
            isString(gffFile, nullOK = TRUE),
            isCharacter(transgeneNames, nullOK = TRUE),
            isCharacter(interestingGroups)
        )
        if (isString(gffFile)) {
            isAFile(gffFile) || isAURL(gffFile)
        }
        h1("bcbioSingleCell")
        alert("Importing bcbio-nextgen single-cell RNA-seq run")
        sampleData <- NULL
        ## Run info ------------------------------------------------------------
        uploadDir <- realpath(uploadDir)
        projectDir <- projectDir(uploadDir)
        sampleDirs <- sampleDirs(uploadDir)
        lanes <- detectLanes(sampleDirs)
        yaml <- import(file.path(projectDir, "project-summary.yaml"))
        dataVersions <-
            importDataVersions(file.path(projectDir, "data_versions.csv"))
        assert(is(dataVersions, "DataFrame"))
        programVersions <-
            importProgramVersions(file.path(projectDir, "programs.txt"))
        assert(is(dataVersions, "DataFrame"))
        log <- import(file.path(projectDir, "bcbio-nextgen.log"))
        ## This step enables our minimal dataset to pass checks.
        tryCatch(
            expr = assert(isCharacter(log)),
            error = function(e) {
                alertWarning(sprintf(
                    "{.file %s} file is empty.",
                    "bcbio-nextgen.log"
                ))
            }
        )
        commandsLog <-
            import(file.path(projectDir, "bcbio-nextgen-commands.log"))
        ## This step enables our minimal dataset to pass checks.
        tryCatch(
            expr = assert(isCharacter(commandsLog)),
            error = function(e) {
                alertWarning(
                    "{.file bcbio-nextgen-commands.log} file is empty."
                )
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
        ## Sample metadata -----------------------------------------------------
        h2("Sample metadata")
        allSamples <- TRUE
        sampleData <- NULL
        if (isString(sampleMetadataFile)) {
            sampleData <- importSampleData(
                file = sampleMetadataFile,
                lanes = lanes,
                pipeline = "bcbio"
            )
            ## Error on incorrect reverse complement input.
            if (isSubset("sequence", colnames(sampleData))) {
                sampleDirSequence <- stri_match_first_regex(
                    str = names(sampleDirs),
                    pattern = "^.+_([ACGT]+)$"
                )[, 2L]
                assert(
                    !identical(
                        sort(sampleDirSequence),
                        sort(as.character(sampleData[["sequence"]]))
                    ),
                    msg = paste(
                        "It appears that the reverse complement sequence of",
                        "the i5 index barcodes were input into the sample",
                        "metadata 'sequence' column. bcbio outputs the revcomp",
                        "into the sample directories, but the forward sequence",
                        "should be used in the R package."
                    )
                )
            }
            ## Allow sample selection by with this file.
            if (nrow(sampleData) < length(sampleDirs)) {
                sampleDirs <- sampleDirs[rownames(sampleData)]
                alert(sprintf(
                    fmt = "Loading a subset of samples: %s.",
                    toInlineString(basename(sampleDirs), n = 5L)
                ))
                allSamples <- FALSE
            }
        }
        ## Assays (counts) -----------------------------------------------------
        h2("Counts")
        ## Note that we're now allowing transcript-level counts.
        counts <- .importCounts(sampleDirs = sampleDirs)
        assert(hasValidDimnames(counts))
        ## Row data (genes/transcripts) ----------------------------------------
        h2("Feature metadata")
        ## Annotation priority:
        ## 1. AnnotationHub.
        ## - Requires `organism` to be declared.
        ## - Ensure that Ensembl release and genome build match.
        ## 2. GTF/GFF file. Use the bcbio GTF if possible.
        ## 3. Fall back to slotting empty ranges. This is offered as support for
        ## complex datasets (e.g. multiple organisms).
        if (isString(organism) && is.numeric(ensemblRelease)) {
            ## AnnotationHub (ensembldb).
            alert("{.fun makeGRangesFromEnsembl}")
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
                alert("{.fun makeGRangesFromGFF}")
                gffFile <- realpath(gffFile)
                rowRanges <- makeGRangesFromGFF(file = gffFile, level = level)
            } else {
                alertWarning("Slotting empty ranges into {.fun rowRanges}.")
                rowRanges <- emptyRanges(rownames(counts))
            }
        }
        assert(is(rowRanges, "GenomicRanges"))
        ## Attempt to get genome build and Ensembl release if not declared.
        ## Note that these will remain NULL when using GTF file (see above).
        if (is.null(genomeBuild)) {
            genomeBuild <- metadata(rowRanges)[["genomeBuild"]]
        }
        if (is.null(ensemblRelease)) {
            ensemblRelease <- metadata(rowRanges)[["ensemblRelease"]]
        }
        ## Column data ---------------------------------------------------------
        h2("Column data")
        colData <- DataFrame(row.names = colnames(counts))
        ## Generate automatic sample metadata, if necessary.
        if (is.null(sampleData)) {
            if (isTRUE(multiplexed)) {
                ## Multiplexed samples without user-defined metadata.
                alertWarning(sprintf(
                    fmt = paste(
                        "{.var %s} is recommended for",
                        "multiplexed samples (e.g. {.val %s})."
                    ),
                    "sampleMetadataFile", umiType
                ))
                sampleData <- minimalSampleData(basename(sampleDirs))
            } else {
                sampleData <- getSampleDataFromYAML(yaml)
            }
        }
        assert(isSubset(rownames(sampleData), names(sampleDirs)))
        ## Join `sampleData` into cell-level `colData`.
        if (identical(nrow(sampleData), 1L)) {
            colData[["sampleId"]] <- as.factor(rownames(sampleData))
        } else {
            colData[["sampleId"]] <- mapCellsToSamples(
                cells = rownames(colData),
                samples = rownames(sampleData)
            )
        }
        sampleData[["sampleId"]] <- as.factor(rownames(sampleData))
        ## Need to ensure the `sampleId` factor levels match up, otherwise we'll
        ## get a warning during the `leftJoin()` call below.
        assert(areSetEqual(
            x = levels(colData[["sampleId"]]),
            y = levels(sampleData[["sampleId"]])
        ))
        levels(sampleData[["sampleId"]]) <- levels(colData[["sampleId"]])
        colData <- leftJoin(colData, sampleData, by = "sampleId")
        assert(
            is(colData, "DataFrame"),
            hasRownames(colData)
        )
        ## Metadata ------------------------------------------------------------
        h2("Metadata")
        cbList <- .importReads(sampleDirs = sampleDirs)
        runDate <- runDate(projectDir)
        interestingGroups <- camelCase(interestingGroups, strict = TRUE)
        assert(isSubset(interestingGroups, colnames(sampleData)))
        metadata <- list(
            "allSamples" = allSamples,
            "bcbioCommandsLog" = commandsLog,
            "bcbioLog" = log,
            "call" = standardizeCall(),
            "cellularBarcodeCutoff" = cutoff,
            "cellularBarcodes" = cbList,
            "dataVersions" = dataVersions,
            "ensemblRelease" = as.integer(ensemblRelease),
            "genomeBuild" = as.character(genomeBuild),
            "gffFile" = as.character(gffFile),
            "interestingGroups" = interestingGroups,
            "lanes" = lanes,
            "level" = level,
            "organism" = as.character(organism),
            "pipeline" = "bcbio",
            "programVersions" = programVersions,
            "projectDir" = projectDir,
            "runDate" = runDate,
            "sampleDirs" = sampleDirs,
            "sampleMetadataFile" = as.character(sampleMetadataFile),
            "umiType" = umiType,
            "uploadDir" = uploadDir,
            "version" = .version,
            "yaml" = yaml
        )
        ## SingleCellExperiment ------------------------------------------------
        object <- makeSingleCellExperiment(
            assays = SimpleList("counts" = counts),
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata,
            transgeneNames = transgeneNames
        )
        ## Return --------------------------------------------------------------
        ## Always prefilter, removing very low quality cells and/or genes.
        object <- calculateMetrics(object = object, prefilter = TRUE)
        ## Bind the `nRead` column into the cell metrics. These are the number
        ## of raw read counts prior to UMI disambiguation that bcbio uses for
        ## initial filtering (`minimum_barcode_depth`` in YAML).
        colData <- colData(object)
        nRead <- .nRead(cbList)
        assert(
            is.integer(nRead),
            isSubset(rownames(colData), names(nRead)),
            areDisjointSets("nRead", colnames(colData))
        )
        colData[["nRead"]] <- unname(nRead[rownames(colData)])
        colData <- colData[, sort(colnames(colData)), drop = FALSE]
        colData(object) <- colData
        bcb <- new(Class = "bcbioSingleCell", object)
        alertSuccess("bcbio single-cell RNA-seq run imported successfully.")
        bcb
    }
