#' Read Sample Metadata File
#'
#' Appends reverse complement sequences by matching against the `sequence`
#' column.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param file Sample metadata file.
#' @param sampleDirs Sample directory paths.
#' @param pipeline Pipeline used to process the samples.
#' @param uniqueNames Ensure sample names are unique.
#'
#' @return [DataFrame].
.readSampleMetadataFile <- function(
    file,
    sampleDirs,
    pipeline) {
    meta <- readFileByExtension(file)

    if (pipeline == "bcbio") {
        # Rename legacy `samplename` column, if set
        if ("samplename" %in% colnames(meta)) {
            warning("Renamed metadata column 'samplename' to 'fileName'")
            meta <- rename(meta, fileName = .data[["samplename"]])
        }

        # Rename `description` to `sampleName`, if set
        if ("description" %in% colnames(meta)) {
            warning("Renamed metadata column 'description' to 'sampleName'")
            meta <- rename(meta, sampleName = .data[["description"]])
        }

        # Check for general required columns
        if (!all(c("fileName", "sampleName") %in% colnames(meta))) {
            stop("'fileName' and 'sampleName' are required", call. = FALSE)
        }

        # Remove incomplete rows
        meta <- meta %>%
            tidy_filter(!is.na(.data[["fileName"]])) %>%
            tidy_filter(!is.na(.data[["sampleName"]]))


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
                                        character(1L)),
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
    tidy_select(meta, metaPriorityCols, everything()) %>%
        as.data.frame %>%
        set_rownames(.[["sampleID"]]) %>%
        as("DataFrame")
}
