# Cell-Cycle Markers
#
# Last updated: 2017-09-09
# Gene annotations: Ensembl Genes 90
#
# This code is derived from:
#   - Tirosh et al, 2015
#   - http://satijalab.org/seurat/cell_cycle_vignette.html

library(devtools)
library(googlesheets)
library(tidyverse)

load_all()

# Allow tidyverse to access Google Sheets
# gs_ls()

# Download the Google sheet (gs)
gs <- gs_key("1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw")

# Get a list of the worksheets (ws)
ws <- gs_ws_ls(gs)
print(ws)

cellCycleMarkers <- lapply(seq_along(ws), function(a) {
    # Stop instead of warn on gene identifier match failure
    tryCatch(
        gs %>%
            gs_read(ws = ws[[a]]) %>%
            mutate(ensgene = symbol2gene(symbol, organism = ws[[a]])) %>%
            select(phase, symbol, ensgene) %>%
            group_by(phase) %>%
            arrange(symbol, .by_group = TRUE),
        warning = function(w) {
            stop(w)
        })
}) %>%
    set_names(ws)

use_data(cellCycleMarkers, compress = "xz", overwrite = TRUE)
