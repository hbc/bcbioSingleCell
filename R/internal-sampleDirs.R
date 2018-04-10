#' Detect Sample Directories
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump makeNames
#' @importFrom bcbioBase sampleDirs
#'
#' @param uploadDir Path to final upload directory.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will abort if no complete sample directories match.
.sampleDirs <- function(
    uploadDir,
    pipeline = c("bcbio", "cellranger")
) {
    assert_all_are_dirs(uploadDir)
    pipeline <- match.arg(pipeline)

    if (pipeline == "bcbio") {
        sampleDirs <- sampleDirs(uploadDir)
    } else if (pipeline == "cellranger") {
        message(paste(
            "CellRanger output directory structure:",
            file.path(
                "<uploadDir>",
                "<sampleName>",
                "outs",
                "filtered_gene_bc_matrices*",
                "<genomeBuild>",
                "matrix.mtx"
            ),
            sep = "\n"
        ))

        matrixFiles <- list.files(
            path = uploadDir,
            pattern = "matrix.mtx",
            include.dirs = FALSE,
            full.names = TRUE,
            recursive = TRUE
        )
        # Subset to only include `filtered_gene_bc_matrices*`. Note that
        # aggregation output is labeled `filtered_gene_bc_matrices_mex` by
        # default.
        matrixFiles <- matrixFiles[grepl(
            x = matrixFiles,
            pattern = file.path(
                paste0("^", uploadDir),
                "[^/]+", # sampleName
                "outs",
                "filtered_gene_bc_matrices([^/]+)?",
                "[^/]+", # genomeBuild
                paste0("matrix.mtx", "$")
            )
        )]

        # Check to ensure that matrices match standardized cellranger export
        if (!length(matrixFiles)) {
            stop("Failed to detect any sample directories")
        }

        # Sample directories nest the matrix files 4 levels deep
        sampleDirs <- matrixFiles %>%
            dirname() %>%
            dirname() %>%
            dirname() %>%
            dirname()
        names(sampleDirs) <- makeNames(basename(sampleDirs), unique = TRUE)

        message(paste(length(sampleDirs), "samples detected"))
    }

    sampleDirs
}
