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
    # Check that uploadDir exists
    if (!dir.exists(uploadDir)) {
        stop("'uploadDir' does not exist")
    }
    uploadDir <- normalizePath(uploadDir)
    if (pipeline == "bcbio") {
        sampleDirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)

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
            make.names(unique = TRUE)
    } else if (pipeline == "cellranger") {
        message(paste(
            "CellRanger output directory structure:",
            file.path("<uploadDir>",
                      "<sampleName>",
                      "outs",
                      "filtered_gene_bc_matrices*",
                      "<genomeBuild>",
                      "matrix.mtx"),
            sep = "\n"
        ))

        matrixFiles <- list.files(
            path = uploadDir,
            pattern = "matrix.mtx",
            include.dirs = FALSE,
            full.names = TRUE,
            recursive = TRUE)
        # Subset to only include `filtered_gene_bc_matrices*`. Note that
        # aggregation output is labeled `filtered_gene_bc_matrices_mex` by
        # default.
        grepl <- grepl(
            x = matrixFiles,
            pattern = file.path(
                paste0("^", uploadDir),
                "[^/]+", # sampleName
                "outs",
                "filtered_gene_bc_matrices[^/]+?",
                "[^/]+", # genomeBuild
                paste0("matrix.mtx", "$")
            )
        )
        matrixFiles <- matrixFiles[grepl]

        # Check to ensure that matrices match standardized cellranger export
        if (!any(matrixFiles)) {
            stop("Failed to match any filtered count matrices", call. = FALSE)
        }

        # Sample directories nest the matrix files 4 levels deep
        sampleDirs <- matrixFiles %>%
            dirname() %>%
            dirname() %>%
            dirname() %>%
            dirname()
        names(sampleDirs) <- make.names(basename(sampleDirs), unique = TRUE)
    } else {
        stop("Unsupported pipeline", call. = FALSE)
    }

    message(paste(length(sampleDirs), "samples detected"))
    sampleDirs
}
