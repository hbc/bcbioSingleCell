library(basejump)
library(pbapply)
library(tidyverse)

deDir <- file.path("results", "differential_expression", "zinger_edger")
dir.create(deDir, recursive = TRUE, showWarnings = FALSE)

# Current Ensembl mouse annotations
anno <- annotable("Mus musculus")

files <-
    dir("data",
        pattern = "zinger_edger_lrt_cluster_",
        full.names = TRUE)

pblapply(seq_along(files), function(a) {
    file <- files[[a]]
    fileStem <- str_replace(basename(file), "\\.[^\\.]+$", "")
    message(fileStem)
    loadAsName(c(lrt = file))
    tbl <- lrt %>%
        # DGELRT object slots the data.frame as `table`
        .[["table"]] %>%
        rownames_to_column("ensgene") %>%
        left_join(anno, by = "ensgene")
    write_csv(tbl, file.path(deDir, paste0(fileStem, ".csv.gz")))
}) %>%
    invisible
