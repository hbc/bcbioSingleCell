#' Load Cell Ranger output
#'
#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) cell
#' counts from `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Path to final upload directory.
#' @param sample_metadata_file (**Required**). Sample barcode metadata file.
#' @param interesting_groups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#'
#' @return [bcbioSCDataSet].
#' @export
load_cellranger <- function(
    upload_dir,
    sample_metadata_file,
    interesting_groups) {
    umi_type <- "chromium"

    # Cell Ranger output files
    barcodes_file <- "barcodes.tsv"
    genes_file <- "genes.tsv"
    matrix_file <- "matrix.mtx"

    # Find sample directories by nested MatrixMarket file
    sample_dirs <- .sample_dirs(upload_dir, nested_file = matrix_file)

    # Check file dependencies
    if (!all(file.exists(file.path(sample_dirs, barcodes_file)))) {
        stop("Barcodes file missing")
    }
    if (!all(file.exists(file.path(sample_dirs, genes_file)))) {
        stop("Genes file missing")
    }
    if (!all(file.exists(file.path(sample_dirs, matrix_file)))) {
        stop("MatrixMarket file missing")
    }

    # Sample metadata
    sample_metadata_file <- normalizePath(sample_metadata_file)
    sample_metadata <- .read_file(sample_metadata_file)

    # Load the Cell Ranger samples
    message("Reading 10X Cell Ranger counts")
    sparse_list <- pblapply(seq_along(sample_dirs), function(a) {
        .sparse_counts(sample_dirs[a])
    }
    ) %>% set_names(names(sample_dirs))
    message("Combining counts into a single sparse matrix")
    sparse_counts <- do.call(cBind, sparse_list)
    rm(sparse_list)

    # Assign genome build based on the Ensembl identifiers. Currently works
    # for human and mouse samples. Improve this method in the future.
    id <- rownames(sparse_counts)[[1L]]
    if (str_detect(id, "^ENSG")) {
        # H. sapiens
        genome_build <- "grch38"
    } else if (str_detect(id, "^ENSMUSG")) {
        # M. musculus
        genome_build <- "grcm38"
        # FIXME Seems like there's a genome build mismatch between Cell Ranger
        # and annotables. For example, `ENSMUSG00000109048`.
    } else {
        stop("Unsupported organism. Need to update package.")
    }


    # Row data ====
    annotable <- .annotable(genome_build)


    # Column data ====
    metrics <- .metrics(sparse_counts, annotable)


    # Metadata ====
    # TODO Consolidate common metadata between bcbio and cellranger
    metadata <- SimpleList(
        upload_dir = upload_dir,
        sample_dirs = sample_dirs,
        interesting_groups = interesting_groups,
        sample_metadata_file = sample_metadata_file,
        sample_metadata = sample_metadata,
        load_date = Sys.Date(),
        umi_type = umi_type,
        wd = getwd(),
        hpc = detect_hpc(),
        session_info = sessionInfo())


    # Package into SummarizedExperiment ====
    se <- .summarized_experiment(
        sparse_counts,
        col_data = metrics,
        row_data = annotable,
        metadata = metadata)
    bcb <- new("bcbioSCDataSet", se)
    bcb
}
