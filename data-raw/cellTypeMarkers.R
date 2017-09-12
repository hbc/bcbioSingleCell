# Cell-Type Markers
#
# Last updated: 2017-09-09
# Gene annotations: Ensembl Genes 90

library(devtools)
library(googlesheets)
library(stringr)
library(tidyverse)

load_all()

# Allow tidyverse to access Google Sheets
# gs_ls()

# Download the Google sheet (gs)
gs <- gs_key("1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0")

# Get a list of the worksheets (ws)
ws <- gs_ws_ls(gs) %>%
    # Remove internal worksheets prefixed with "_"
    .[!str_detect(., "^_")]
print(ws)

cellTypeMarkers <- lapply(seq_along(ws), function(a) {
    # Stop instead of warn on gene identifier match failure
    tryCatch(
        gs %>%
            gs_read(ws = ws[[a]]) %>%
            group_by(symbol) %>%
            nest %>%
            mutate(ensgene = symbol2gene(symbol, organism = ws[[a]])) %>%
            unnest %>%
            select(cell, symbol, ensgene) %>%
            group_by(cell) %>%
            arrange(symbol, .by_group = TRUE),
        warning = function(w) {
            stop(w)
        })
}) %>%
    set_names(ws)

use_data(cellTypeMarkers, compress = "xz", overwrite = TRUE)
