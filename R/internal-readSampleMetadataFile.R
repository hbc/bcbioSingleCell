#' Read Sample Metadata File
#'
#' Appends reverse complement sequences by matching against the `sequence`
#' column.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param file Sample metadata file.
#' @param sampleDirs Sample directory paths.
#' @param pipeline Pipeline used to process the samples.
#' @param uniqueNames Ensure sample names are unique.
#'
#' @return [data.frame].
.readSampleMetadataFile <- function(
    file,
    sampleDirs,
    pipeline) {
    meta <- readFileByExtension(file)

    if (pipeline == "bcbio") {
        # Stop on legacy `samplename` column. We need to work on improving the
        # consistency in examples or the internal handlng of file and sample
        # names in a future update.
        if ("samplename" %in% colnames(meta)) {
            stop(paste(
                "'samplename' is used in some bcbio examples to define FASTQ",
                "file names, and 'description' to define sample names. Here",
                "we are using 'sampleName' (note lowerCamelCase) in the",
                "package to define sample names. When passing in a sample",
                "metadata file here, use 'fileName' for file names",
                "(e.g. 'control_replicate_1.fastq.gz') and either",
                "'description' or 'sampleName' for sample names",
                "(e.g. 'control replicate 1')"
            ))
        }

        # Rename `description` to `sampleName`, if `sampleName` is unset
        if ("description" %in% colnames(meta) &
            !"sampleName" %in% colnames(meta)) {
            message("Renamed metadata column 'description' to 'sampleName'",
                    call. = FALSE)
            meta <- dplyr::rename(meta, sampleName = .data[["description"]])
        }

        # Check for general required columns
        if (!all(c("fileName", "sampleName") %in% colnames(meta))) {
            stop("'fileName' and 'sampleName' are required", call. = FALSE)
        }

        # Remove incomplete rows
        meta <- meta %>%
            dplyr::filter(!is.na(.data[["fileName"]])) %>%
            dplyr::filter(!is.na(.data[["sampleName"]]))

        # Check that sample names are unique
        if (any(duplicated(meta[["sampleName"]]))) {
            stop("'sampleName' column does not contain unique values",
                 call. = FALSE)
        }

        # Check if samples are demultiplexed
        if (length(unique(meta[["fileName"]])) == nrow(meta)) {
            # SureCell
            demultiplexed <- TRUE
        } else {
            # inDrop
            demultiplexed <- FALSE
        }

        if (isTRUE(demultiplexed)) {
            meta[["sampleID"]] <- meta[["sampleName"]]
        } else {
            # Match the UMI demultiplexed sample directories (e.g. inDrop)
            if (!"sequence" %in% colnames(meta)) {
                stop("Index i7 sequence required to generate 'sampleID'")
            }
            meta <- meta %>%
                mutate(revcomp = vapply(.data[["sequence"]],
                                        revcomp,
                                        character(1)),
                       sampleID = paste(
                           make.names(.data[["fileName"]]),
                           .data[["revcomp"]],
                           sep = "_"))
        }
    } else if (pipeline == "cellranger") {
        # Check for general required columns
        if (!all(c("fileName", "sampleName") %in% colnames(meta))) {
            stop("'fileName' and 'sampleName' are required", call. = FALSE)
        }
        meta[["sampleID"]] <- make.names(meta[["fileName"]])
    } else {
        stop("Unsupported pipeline")
    }

    # Arrange rows by `sampleID`
    meta <- arrange(meta, !!sym("sampleID"))

    # Check that `sampleID` matches `sampleDirs`
    if (!all(meta[["sampleID"]] %in% names(sampleDirs))) {
        stop("Sample directory names don't match the sample metadata file")
    }

    # Return
    dplyr::select(meta, metaPriorityCols, everything()) %>%
        as.data.frame() %>%
        set_rownames(.[["sampleID"]])
}
