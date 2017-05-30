read_bcbio_barcodes <- function(run) {
    files <- run$sample_dirs %>%
        file.path(., paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        set_names(names(run$sample_dirs))
    if (!all(file.exists(files))) {
        stop("Barcode TSV file missing")
    }
    list <- mclapply(seq_along(files), function(a) {
        read_barcode_file(files[a])
    }) %>% set_names(names(files))
}
