#' Detect Sample Directories
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param uploadDir Upload directory.
#' @param pipeline Pipeline used to generate the samples.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [stop()] if no complete sample directories match.
.sampleDirs <- function(uploadDir, pipeline = "bcbio") {
    if (pipeline == "bcbio") {
        sampleDirs <- list.dirs(
            uploadDir, full.names = TRUE, recursive = FALSE)
        # Remove the nested `projectDir`
        if (any(str_detect(basename(sampleDirs), projectDirPattern))) {
            sampleDirs <- sampleDirs %>%
                .[!str_detect(basename(.), projectDirPattern)]
        }
        if (length(sampleDirs) == 0L) {
            stop("Failed to detect any sample directories",
                 call. = FALSE)
        }
        # Sample directory: `sampleName`-`barcode`
        split <- basename(sampleDirs) %>%
            str_split("-")
        sampleNames <-
            sapply(seq_along(split), function(a) {
            paste(camel(split[[a]][[1L]]),
                  split[[a]][[2L]],
                  sep = "_")
        })
        names(sampleDirs) <- sampleNames
    } else if (pipeline == "cellranger") {
        # Faster directory matching
        subdirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)
        matches <- dir.exists(file.path(
            subdirs, "outs", "filtered_gene_bc_matrices"))
        if (!any(matches)) {
            stop("Failed to detect 'filtered_gene_bc_matrices'",
                 call. = FALSE)
        }
        sampleDirs <- subdirs[matches]
        names(sampleDirs) <- camel(basename(sampleDirs))
    } else {
        stop("Unsupported pipeline")
    }

    message(paste(length(sampleDirs), "samples detected"))
    sampleDirs
}
