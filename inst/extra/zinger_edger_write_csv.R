de_dir <- file.path("results", "differential_expression", "zinger_edger")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

# Current Ensembl mouse annotations
anno <- basejump::annotable("Mus musculus")

cluster_de_files <-
    dir("data",
        pattern = "zinger_edger_lrt_cluster_",
        full.names = TRUE)

pbapply::pblapply(seq_along(cluster_de_files), function(a) {
    dge_name <- load(cluster_de_files[[a]])
    message(dge_name)
    dge <- get(dge_name, inherits = FALSE)
    tbl <- dge$table %>%
        tibble::rownames_to_column("ensgene") %>%
        dplyr::left_join(anno, by = "ensgene")
    readr::write_csv(tbl, file.path(de_dir, paste0(dge_name, ".csv.gz")))
}) %>%
    invisible
