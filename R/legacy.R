#' Create a skeleton bcbio project
#'
#' @keywords internal
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param base_dir Base directory in which to create skeleton bcbio structure.
#' @param config_dir Configuration directory.
#' @param upload_dir Upload directory.
#' @param project_dir Project directory (YYYY-MM-DD_template), nested inside
#'   upload directory.
#' @param sample_dirs Sample directories, nested inside upload directory.
#' @param organism Organism to use.
#'
#' @return bcbio run skeleton.
#' @export
create_bcbio_skeleton <- function(
    base_dir = getwd(),
    run_dir = "skeleton",
    config_dir = "config",
    upload_dir = "final",
    project_dir = NULL,
    sample_dirs = NULL,
    organism,
    intgroup = "sample_name") {
    # run_dir
    if (!is.null(run_dir)) {
        run_dir <- file.path(base_dir, run_dir) %>% normalizePath
        dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    } else {
        stop("Run directory required")
    }

    # config_dir
    if (!is.null(config_dir)) {
        config_dir <- file.path(run_dir, config_dir)
        dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
    } else {
        stop("Configuration directory required")
    }

    # upload_dir
    if (!is.null(upload_dir)) {
        upload_dir <- file.path(run_dir, upload_dir)
        dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)
    } else {
        stop("Final upload directory required")
    }

    # project_dir
    if (!is.null(project_dir)) {
        project_dir <- file.path(upload_dir, project_dir)
        dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # sample_dirs
    if (!is.null(sample_dirs)) {
        sample_dirs <- file.path(upload_dir, sample_dirs)
        dir.create(sample_dirs, recursive = TRUE, showWarnings = FALSE)
    }

    list(run_dir = run_dir,
         config_dir = config_dir,
         upload_dir = upload_dir,
         project_dir = project_dir,
         sample_dirs = sample_dirs,
         organism = organism,
         intgroup = intgroup)
}



# Read a MatrixMart file, setting the row and column names
#
# Note that for a bcbio run, this function will return transcript-level counts.
#
# @author Rory Kirchner
# @author Michael Steinbaugh
#
# @param matrix_file MatrixMart file to read.
#
# @return Sparse counts matrix.
.counts2 <- function(matrix_file) {
    # Detect gzip file
    pattern <- "\\.gz$"
    if (str_detect(matrix_file, pattern)) {
        matrix_file <- gunzip(matrix_file)
    }

    row_file <- paste0(matrix_file, ".rownames")
    col_file <- paste0(matrix_file, ".colnames")

    if (!file.exists(matrix_file)) {
        stop("MatrixMart counts file missing")
    }
    if (!file.exists(row_file)) {
        stop("Rownames file could not be found")
    }
    if (!file.exists(col_file)) {
        stop("Colnames file could not be found")
    }

    counts <- readMM(matrix_file)
    rownames(counts) <- read_lines(row_file)
    colnames(counts) <- read_lines(col_file)

    # [fix] Correct malformed celluar barcodes.
    # Need to update on the bcbio-nextgen platform side.
    if (any(str_detect(colnames(counts), "\\:[ACGT]{16}$"))) {
        colnames(counts) <- colnames(counts) %>%
            str_replace("\\:([ACGT]{8})([ACGT]{8})$", "\\:\\1-\\2")
    }

    # Coerce dgTMatrix to dgCMatrix
    as(counts, "dgCMatrix")
}
