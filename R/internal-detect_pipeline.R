#' Detect scRNA-seq analysis pipeline
#'
#' Currently supports bcbio-nextgen and 10X Chromium.
#'
#' @rdname detect_pipeline
#' @author Michael Steinbaugh
#'
#' @param upload_dir Final upload directory.
#'
#' @return logical.
.detect_pipeline <- function(upload_dir) {
    matrices <- list.files(
        upload_dir, pattern = "*.mtx$",
        full.names = TRUE, recursive = TRUE) %>%
        normalizePath %>%
        set_names(basename(.))
    if (!length(matrices)) {
        stop("No count matrices detected")
    }
    parent_dirs <- dirname(matrices)
    if (any(str_detect(names(matrices), "^tagcounts\\.mtx$"))) {
        # bcbio-nextgen ====
        # Check for barcode (`.colnames`) and transcript (`.rownames`)
        # dependency files
        if (all(file.exists(paste0(matrices, ".colnames"))) &
            all(file.exists(paste0(matrices, ".rownames")))) {
            "bcbio"
        } else {
            NULL
        }
    } else if (all(str_detect(names(matrices), "^matrix\\.mtx$"))) {
        # 10X Chromium CellRanger ====
        # Check for barcode (`barcodes.tsv`) and gene (`genes.tsv`) dependency
        # files
        if (all(file.exists(file.path(parent_dirs, "barcodes.tsv"))) &
            all(file.exists(file.path(parent_dirs, "genes.tsv")))) {
            "cellranger"
        } else {
            NULL
        }
    } else {
        NULL
    }
}
