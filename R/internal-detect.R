#' Detect scRNA-seq analysis pipeline
#'
#' Currently supports bcbio-nextgen and 10X Genomics Cell Ranger.
#'
#' @rdname detect_pipeline
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param sample_dirs Sample directories.
#'
#' @return Pipeline string. Stops on detection failure.
.detect_pipeline <- function(sample_dirs) {
    bcbio_matrix <-
        file.path(sample_dirs,
                  str_c(basename(sample_dirs), ".mtx"))
    cellranger_matrix <-
        file.path(sample_dirs, "matrix.mtx")
    if (all(file.exists(bcbio_matrix))) {
        # bcbio-nextgen ====
        # Check for barcode (`.colnames`) and transcript (`.rownames`)
        # dependency files
        if (all(file.exists(str_c(bcbio_matrix, ".colnames"))) &
            all(file.exists(str_c(bcbio_matrix, ".rownames")))) {
            pipeline <- "bcbio-nextgen"
        } else {
            pipeline <- NULL
        }
    } else if (all(file.exists(cellranger_matrix))) {
        # 10X Chromium CellRanger ====
        # Check for barcode (`barcodes.tsv`) and gene (`genes.tsv`) dependency
        # files
        if (all(file.exists(file.path(sample_dirs, "barcodes.tsv"))) &
            all(file.exists(file.path(sample_dirs, "genes.tsv")))) {
            pipeline <- "cellranger"
        } else {
            pipeline <- NULL
        }
    } else {
        pipeline <- NULL
    }
    if (is.null(pipeline)) {
        stop("Pipeline detection failed")
    }
    pipeline
}
