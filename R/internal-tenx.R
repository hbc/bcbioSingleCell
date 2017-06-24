#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) output
#'
#' Read into sparse matrix from barcodes.tsv, genes.tsv, and matrix.mtx files.
#'
#' @rdname tenx
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Final upload directory.
#' @param aggregate Aggregate samples into single sparse matrix.
#'
#' @return Sparse counts matrix (`dgCMatrix`).
.cellranger <- function(upload_dir, aggregate = TRUE) {
    # Recurse through specific data directory and identify sample subdirectories
    # by the presence of a `matrix.mtx` counts file.
    sample_dirs <- list.files(
        upload_dir,
        full.names = TRUE,
        pattern = "matrix.mtx",
        recursive = TRUE) %>%
        normalizePath %>%
        dirname %>%
        set_names(basename(.))
    message(paste(length(sample_dirs), "samples detected"))
    message("Reading 10X Chromium samples...")
    lst <- pblapply(seq_along(sample_dirs), function(a) {
        name <- names(sample_dirs)[a]
        sample_dir <- sample_dirs[a]
        barcodes_file <- file.path(sample_dir, "barcodes.tsv")
        if (!file.exists(barcodes_file)) {
            stop("Barcodes file missing")
        }
        genes_file <- file.path(sample_dir, "genes.tsv")
        if (!file.exists(genes_file)) {
            stop("Gene names file missing")
        }
        matrix_file <- file.path(sample_dir, "matrix.mtx")
        if (!file.exists(matrix_file)) {
            stop("MatrixMart counts file missing")
        }

        # Read MatrixMart expression counts matrix
        counts <- readMM(matrix_file)

        # Assign barcodes to colnames
        barcodes <- read_lines(barcodes_file)
        # Check for `-1` suffix in barcode names, indicative of 10X data
        pattern <- "\\-1$"
        if (all(grepl(pattern, barcodes))) {
            barcodes <- str_replace(barcodes, pattern, "")
        }
        # Append the sample name to the barcodes. This will be used after
        # list is generated to [cBind()] all counts to a single sparse matrix.
        colnames(counts) <- paste(name, barcodes, sep = "-")

        # Assign gene names (symbols) to rownames
        genes <- read_tsv(genes_file,
                          col_names = c("ensembl_gene_id",
                                        "external_gene_name"),
                          col_types = cols())
        rownames(counts) <- make.unique(genes$external_gene_name)

        # Coerce dgTMatrix to dgCMatrix (sparse)
        as(counts, "dgCMatrix")
    }) %>% set_names(names(sample_dirs))

    # TODO Offload this step to a more general internal aggregation function?
    if (isTRUE(aggregate)) {
        # Aggregate the counts, if desired
        message("Aggregating samples to a single sparse matrix...")
        do.call(cBind, lst)
    } else {
        # Otherwise, return a per-sample list
        lst
    }
}
