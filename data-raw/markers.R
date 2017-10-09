# Cell Markers
#
# Last updated: 2017-09-13
# Gene annotations: Ensembl Genes 90
#
# This code is derived from:
#   - Tirosh et al, 2015
#   - http://satijalab.org/seurat/cell_cycle_vignette.html

library(basejump)
library(devtools)
library(googlesheets)
library(stringr)
library(tidyverse)

load_all()

# Ensembl release version
# bioc-devel required for release 90 annotations (9/13/17)
# Use `BiocInstaller::useDevel()` and then run `biocLite()`
release <- 90

# Here we're matching the stored Ensembl identifiers (`ensgene`) using
# ensembldb to obtain the latest symbol names from Ensembl.

# Allow tidyverse to access Google Sheets
gs_ls()

# Cell Cycle Markers ===========================================================
# Download the Google sheet (gs)
gs <- gs_key("1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw")

# Get a list of the worksheets (ws)
ws <- gs_ws_ls(gs)
print(ws)

cellCycleMarkers <- lapply(seq_along(ws), function(a) {
    gs %>%
        gs_read(ws = ws[[a]]) %>%
        dplyr::select(phase, ensgene) %>%
        mutate(symbol = gene2symbol(
            ensgene, organism = ws[[a]], release = release)) %>%
        group_by(phase) %>%
        arrange(symbol, .by_group = TRUE)
}) %>%
    set_names(ws)

# Cell Type Markers ============================================================
# Download the Google sheet (gs)
gs <- gs_key("1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0")

# Get a list of the worksheets (ws)
ws <- gs_ws_ls(gs) %>%
    # Remove internal worksheets prefixed with "_"
    .[!str_detect(., "^_")]
print(ws)

cellTypeMarkers <- lapply(seq_along(ws), function(a) {
    gs %>%
        gs_read(ws = ws[[a]]) %>%
        dplyr::select(cell, ensgene) %>%
        mutate(symbol = gene2symbol(
            ensgene, organism = ws[[a]], release = release)) %>%
        group_by(cell) %>%
        arrange(symbol, .by_group = TRUE)
}) %>%
    set_names(ws)

# Save RData ===================================================================
use_data(
    cellCycleMarkers,
    cellTypeMarkers,
    compress = "xz",
    overwrite = TRUE)
