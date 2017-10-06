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
#'
#' @return [data.frame].
.readSampleMetadataFile <- function(file) {
    meta <- readFileByExtension(file)

    # Stop on legacy `samplename` column. We need to work on improving the
    # consistency in examples or the internal handlng of file and sample
    # names in a future update.
    if ("samplename" %in% colnames(meta)) {
        stop(paste(
            "'samplename' is used in some bcbio examples for FASTQ file names,",
            "and 'description' for sample names. Here we are using 'fileName'",
            "for FASTQ file names (e.g. 'control_replicate_1.fastq.gz'),",
            "'description' for multiplexed per file sample names",
            "(e.g. 'control replicate 1', and 'sampleName' for multiplexed",
            "sample names (i.e. inDrop barcoded samples)."
        ))
    }

    if (any(duplicated(meta[["description"]]))) {
        multiplexedFASTQ <- TRUE
    } else {
        multiplexedFASTQ <- FALSE
    }

    if (isTRUE(multiplexedFASTQ)) {
        requiredCols <- c("fileName", "description", "sampleName", "sequence")
        if (!all(requiredCols %in% colnames(meta))) {
            stop(paste("Required columns:", toString(requiredCols)))
        }
    } else {
        requiredCols <- c("fileName", "description")
        if (!all(requiredCols %in% colnames(meta))) {
            stop(paste("Required columns:", toString(requiredCols)))
        }
        # Check for duplicate `description` and `sampleName`
        if (all(c("description", "sampleName"))) {
            stop(paste(
                "Specify only 'description' and omit 'sampleName' for",
                "demultiplexed FASTQ file metadata"
            ))
        }
        meta[["sampleName"]] <- meta[["description"]]
    }

    # Remove incomplete rows
    meta <- meta %>%
        dplyr::filter(!is.na(.data[["description"]])) %>%
        dplyr::filter(!is.na(.data[["sampleName"]]))

    # Check that sample names are unique
    if (any(duplicated(meta[["sampleName"]]))) {
        stop("Sample name descriptions are not unique",
             call. = FALSE)
    }

    # Set the `sampleID` column
    if (isTRUE(multiplexedFASTQ)) {
        # The per sample directories are created by combining the
        # `sampleName` column with the reverse complement (`revcomp`) of the
        # index barcode sequence (`sequence`)
        meta <- meta %>%
            mutate(
                revcomp = vapply(.data[["sequence"]], revcomp, character(1)),
                sampleID = paste(
                    .data[["description"]],
                    .data[["revcomp"]],
                    sep = "_")
            )
    } else {
        meta[["sampleID"]] <- meta[["sampleName"]]
    }

    meta %>%
        mutate(
            revcomp = vapply(.data[["sequence"]], revcomp, character(1)),
            # Match the sample directories exactly
            sampleID = paste(
                .data[["description"]],
                .data[["revcomp"]],
                sep = "-"),
            # Now sanitize following the sample rules in `.sampleDirs()`
            sampleID = make.names(
                str_replace_all(.data[["sampleID"]], "-", "_"))
        ) %>%
        dplyr::select(metaPriorityCols, everything()) %>%
        arrange(!!sym("sampleID")) %>%
        as.data.frame() %>%
        set_rownames(.[["sampleID"]])
}
