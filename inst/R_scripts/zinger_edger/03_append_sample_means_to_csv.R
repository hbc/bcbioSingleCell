library(basejump)
library(Matrix)
library(pbapply)
library(stringr)
library(tidyverse)

# Load up the Seurat object
loadData(seurat)

degTableDir <- file.path(
    "results",
    "differential_expression",
    "zinger_edger"
)
dir.exists(degTableDir)

newTableDir <- file.path(
    "results",
    "differential_expression",
    "zinger_edger_with_counts"
)
dir.create(newTableDir, recursive = TRUE, showWarnings = FALSE)

# Let's get a list of the CSV files per cluster
csv <- dir(degTableDir, pattern = "*.csv.gz", full.names = TRUE)

# Now we can loop across the CSV files
pblapply(seq_along(csv), function(a) {
    file <- csv[[a]]
    message(basename(file))
    # Get the cluster ident from the file name
    ident <- str_match(file, "cluster_([0-9]+)\\.csv") %>%
        .[, 2L]

    # Subset the original seurat object by the ident
    # We need to fetch the counts by cluster
    subset <- SubsetData(seurat, ident.use = ident)

    # Calculate the mean expression per cell per sample
    samples <- unique(subset@meta.data[["sampleName"]])

    # Run the nested sample mean loop
    means <- lapply(seq_along(samples), function(b) {
        cells <- subset@meta.data %>%
            .[.[["sampleName"]] == samples[[b]], ] %>%
            rownames()
        # The log normalized counts from Seurat are in the `@data` slot. Don't
        # use `@raw.data` or `@scale.data`.
        counts <- subset@data[, cells, drop = FALSE]
        # Substitute the geneID for symbol in the rownames. We'll use this
        # later for easier joining back to the zingeR/edgeR CSV
        rownames(counts) <- names(rownames(counts))
        # Calculate the mean expression across the cells. If there is only
        # 1 cell present, then just return the values.
        if (ncol(counts) == 1L) {
            as.numeric(counts)
        } else {
            Matrix::rowMeans(counts)
        }
    }) %>%
        set_names(samples) %>%
        as.data.frame() %>%
        rownames_to_column("geneID")

    # Load and subset the data frame. The annotated CSVs we're loading haven't
    # been filtered by P value.
    df <- read_csv(file) %>%
        # Remove rows without an adjusted P value
        dplyr::filter(!is.na(padjFilter)) %>%
        # Remove rows with an adjusted P value >= 0.05
        dplyr::filter(padjFilter < 0.05) %>%
        # Now arrange by adjusted P value
        arrange(padjFilter, geneName) %>%
        left_join(means, by = "geneID")
    write_csv(df, file.path(newTableDir, basename(file)))
}) %>%
    invisible()
