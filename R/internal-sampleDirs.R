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
.sampleDirs <- function(
    uploadDir,
    pipeline = "bcbio") {
    if (pipeline == "bcbio") {
        sampleDirs <- list.dirs(
            uploadDir, full.names = TRUE, recursive = FALSE)
        # Remove the nested `projectDir`
        if (any(grepl(x = basename(sampleDirs),
                      pattern = projectDirPattern))) {
            sampleDirs <- sampleDirs %>%
                .[!grepl(x = basename(.),
                         pattern = projectDirPattern)]
        }
        if (length(sampleDirs) == 0) {
            stop("Failed to detect any sample directories",
                 call. = FALSE)
        }
        names(sampleDirs) <- basename(sampleDirs) %>%
            gsub(x = .,
                 pattern = "-",
                 replacement = "_") %>%
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
