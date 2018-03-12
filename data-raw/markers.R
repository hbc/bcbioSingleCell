# Cell Markers
#
# Last updated: 2018-03-11
# Gene annotations: Ensembl Genes 90
#
# This code is derived from:
#   - Tirosh et al, 2015
#   - http://satijalab.org/seurat/cell_cycle_vignette.html

devtools::load_all()
library(googlesheets)
library(tidyverse)

# Ensembl release version
release <- 90L

# Here we're matching the stored Ensembl identifiers (`geneID`) using
# ensembldb to obtain the latest symbol names from Ensembl.

# Allow tidyverse to access Google Sheets
gs_ls()

# Cell cycle markers ===========================================================
# Download the Google sheet (gs)
gs <- gs_key("1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw")

# Get a list of the worksheets (ws)
ws <- gs_ws_ls(gs)
print(ws)

cellCycleMarkers <- lapply(ws, function(ws) {
    gs %>%
        gs_read(ws = ws) %>%
        select(phase, geneID) %>%
        mutate(
            geneName = convertGenesToSymbols(
                geneID,
                organism = ws,
                release = release)
        ) %>%
        group_by(phase) %>%
        arrange(geneName, .by_group = TRUE)
})
names(cellCycleMarkers) <- camel(ws)

# Cell type markers ============================================================
# Download the Google sheet (gs)
gs <- gs_key("1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0")

# Get a list of the worksheets (ws)
ws <- gs_ws_ls(gs) %>%
    # Remove internal worksheets prefixed with "_"
    .[!str_detect(., "^_")]
print(ws)

cellTypeMarkers <- lapply(ws, function(ws) {
    gs %>%
        gs_read(ws = ws) %>%
        dplyr::select(cell, geneID) %>%
        mutate(
            geneName = convertGenesToSymbols(
                geneID,
                organism = ws[[a]],
                release = release)
        ) %>%
        group_by(cell) %>%
        arrange(geneName, .by_group = TRUE)
})
names(cellTypeMarkers) <- camel(ws)

# Save R data ==================================================================
devtools::use_data(
    cellCycleMarkers,
    cellTypeMarkers,
    compress = "xz",
    overwrite = TRUE
)
