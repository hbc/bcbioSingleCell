#' Detect Sample Directories
#'
#' @author Michael Steinbaugh
#'
#' @param uploadDir Upload directory.
#' @param pipeline Pipeline used to generate the samples.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will [stop()] if no complete sample directories match.
#' @noRd
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
        names(sampleDirs) <- basename(sampleDirs) %>%
            str_replace_all("-", "_") %>%
            make.names()
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
        names(sampleDirs) <- make.names(basename(sampleDirs))
    } else {
        stop("Unsupported pipeline")
    }

    message(paste(length(sampleDirs), "samples detected"))
    sampleDirs
}
