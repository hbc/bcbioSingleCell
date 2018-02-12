# Cell Markers
#
# Last updated: 2017-11-09
# Gene annotations: Ensembl Genes 90
#
# This code is derived from:
#   - Tirosh et al, 2015
#   - http://satijalab.org/seurat/cell_cycle_vignette.html

devtools::load_all()
library(googlesheets)
library(tidyverse)

# Ensembl release version
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
        mutate(
            symbol = convertGenesToSymbols(
                ensgene,
                organism = ws[[a]],
                release = release)
        ) %>%
        group_by(phase) %>%
        arrange(symbol, .by_group = TRUE)
})
names(cellCycleMarkers) <- camel(ws)

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
        mutate(
            symbol = convertGenesToSymbols(
                ensgene,
                organism = ws[[a]],
                release = release)
        ) %>%
        group_by(cell) %>%
        arrange(symbol, .by_group = TRUE)
})
names(cellTypeMarkers) <- camel(ws)

# Save RData ===================================================================
devtools::use_data(
    cellCycleMarkers,
    cellTypeMarkers,
    compress = "xz",
    overwrite = TRUE)
