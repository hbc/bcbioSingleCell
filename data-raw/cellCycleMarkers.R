# Cell-Cycle Markers
# Last compiled on September 7, 2017 using Ensembl 90 (annotables)
# @seealso
#   - Tirosh et al, 2015
#   - http://satijalab.org/seurat/cell_cycle_vignette.html
library(devtools)
library(tidyverse)
library(readxl)

load_all()

# Store all the sheets in a list
file <- file.path("data-raw", "cellCycleMarkers.xlsx")
sheets <- excel_sheets(file)

cellCycleMarkers <- lapply(seq_along(sheets), function(a) {
    organism <- sheets[[a]]
    # Stop instead of warn on gene identifier match failure
    tryCatch(
        read_excel(file, sheet = organism) %>%
            as.data.frame %>%
            set_rownames(.$symbol) %>%
            symbol2gene(organism = organism) %>%
            rownames_to_column("ensgene") %>%
            as_tibble %>%
            select(phase, ensgene, symbol) %>%
            group_by(phase) %>%
            arrange(symbol, .by_group = TRUE),
        warning = function(w) {
            stop(w)
        })
}) %>%
    set_names(sheets)

use_data(cellCycleMarkers, compress = "xz", overwrite = TRUE)
